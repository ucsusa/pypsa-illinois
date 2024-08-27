import numpy as np
import sys
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
from unyt import m, s, MW, W, kg, g

sys.path.append("functions")

from nrel_data_api import parameters, make_csv_url

model_years = np.array(snakemake.config['model_years']).astype('int')

def handle_datetime(dataframe, use_actual=False):
    """
    Combines time columns into a single timestamp column.
    Expects columns ['year','month','day','hour'].

    Parameters
    ----------
    dataframe : :class:`pd.DataFrame`
        A pandas dataframe.
    """
    frame = dataframe.copy()
    
    if use_actual:
            cols=['Year','Month', 'Day', 'Hour', 'Minute']
            frame = frame.T.drop_duplicates().T
            frame['time'] = pd.to_datetime(frame[cols])
            frame.drop(columns=cols, inplace=True)
            frame.set_index('time', inplace=True)
    
    else:
        timestamps = pd.DatetimeIndex([])
        for year in model_years:
            period = pd.date_range(start=f"{year}-01-01",
                                            freq=f"1h",
                                            periods=8760)
            timestamps = timestamps.append(period)
        
        frame.index = timestamps
        try:
            frame.set_index(timestamps,inplace=True)
            frame.drop(columns=['Year','Month','Day','Hour','Minute'],inplace=True)
        except:
            raise ValueError
    
    return frame
    

def retrieve_solar_timeseries(region, use_actual=False):
    """
    Retrieves data from NREL's national solar radiation database (NSRDB).

    Parameters
    ----------
    region : :class:`gpd.GeoDataFrame`
        A geopandas dataframe containing modeled bus regions.
    """
   
    parameters['attr_list'] = ['ghi']
    outer_pbar = tqdm(snakemake.config['solar_years'], position=0, leave=True)
    all_frames = []
    for year in outer_pbar:
        outer_pbar.set_description(f"Processing {year}")
        parameters['year'] = int(year)
        frames = []
        inner_pbar = tqdm(region[['name','x','y']].values, position=1, leave=True)
        for n, i, j in inner_pbar:
            inner_pbar.set_description(f"Processing {n}")
            parameters['lon'] = i
            parameters['lat'] = j
            URL = make_csv_url(parameters=parameters, 
                            kind='solar')
            df = pd.read_csv(URL, skiprows=2)[:8760]
            df.rename(columns={'GHI':f"{n}"}, inplace=True)
            frames.append(df)
            
            solar_df = pd.concat(frames, axis=1)
    
        all_frames.append(solar_df)
    full_df = pd.concat(all_frames, axis=0)
    full_df = handle_datetime(full_df, use_actual=use_actual)
    
    return full_df


def retrieve_wind_timeseries(region, use_actual=False):
    """
    Retrieves data from NREL's national solar radiation database (NSRDB).

    Parameters
    ----------
    region : :class:`gpd.GeoDataFrame`
        A geopandas dataframe containing modeled bus regions.
    """
    wind_attr = ['windspeed_80m']
    parameters['attr_list'] = wind_attr
    
    outer_pbar = tqdm(snakemake.config['wind_years'], leave=True, position=0)
    all_frames = []
    for year in outer_pbar:
        outer_pbar.set_description(f"Processing {year}")
        
        parameters['year'] = int(year)
        pbar = tqdm(region[['name','x','y']].values, position=1, leave=True)
        frames = []
        for n, i, j in pbar:
            pbar.set_description(f"Processing {n}")
            parameters['lon'] = i
            parameters['lat'] = j
            URL = make_csv_url(parameters=parameters, 
                            kind='wind')
            df = pd.read_csv(URL, skiprows=1)[:8760]
            df.rename(columns={'wind speed at 80m (m/s)':f"{n}"}, inplace=True)
            df.to_csv(f"data/timeseries/wind/{n}_wind_{year}.csv")
            frames.append(df)
        
        wind_df = pd.concat(frames, axis=1)
        all_frames.append(wind_df)
    full_df = pd.concat(all_frames, axis=0)
    full_df = handle_datetime(full_df, use_actual=use_actual)
    
    return full_df


def turbine_power(v):
    """
    Calculates the power output of a wind turbine.
    
    Parameters
    ----------
    v : float
        Wind speed in meters per second (m/s).
    Returns
    -------
    power : float
        Power in megawatts (MW).
    """
    turbine_params = snakemake.config['turbine_params']
    
    cut_in = float(turbine_params['cut_in'])*m/s
    cut_out = float(turbine_params['cut_out'])*m/s
    rated = float(turbine_params['rated'])*m/s
    diameter = float(turbine_params['diameter'])*m
    rated_power = float(turbine_params['rated_power'])*MW
    air_density = float(turbine_params['air_density'])*kg/m**3
    
    
    power = lambda v: np.min((((0.5*np.pi/4)*(air_density*diameter**2)*v**3).to(MW), 
                              rated_power))*MW
    
    wind_speed = v*m/s
    
    if wind_speed < cut_in:
        return 0*MW
    elif (wind_speed > rated) and (wind_speed < cut_out):
        return rated_power
    elif wind_speed >= cut_out:
        return 0*MW
    elif wind_speed >= cut_in:
        return power(wind_speed)


def process_wind_timeseries(df, normalize=True):
    """
    Converts windspeed (m/s) timeseries data to a 
    hypothetical power production given some wind turbine
    parameters (specified in `config.yml`).

    Parameters
    ----------
    df : :class:pd.DataFrame
        A pandas dataframe with a column of windspeed data
        in (m/s).
    normalize : boolean
        Whether the data should be normalized. Default is true.
    """
    frame = df.copy()
    
    for col in frame.columns:
        frame[col] = frame[col].apply(turbine_power)
        
    
    if normalize:
        frame = frame.divide(frame.max(axis=0), axis=1)
    
    return frame 
    
    
def process_solar_timeseries(df, normalize=True):
    """
    Converts solar radiation timeseries to a
    hypothetical power production.

    Parameters
    ----------
    df : _type_
        _description_
    normalize : bool, optional
        Whether the data should be normalized. Default is true.
    """
    frame = df.copy()
    if normalize:
            frame = frame.divide(frame.max(axis=0), axis=1)
        
    return frame 

if __name__ == "__main__":
    
    regions = gpd.read_file(snakemake.input.supply_regions)
    use_actual = snakemake.config['use_actual']
    
    # solar data
    df = retrieve_solar_timeseries(regions, use_actual=use_actual)
    df = process_solar_timeseries(df)
    df.to_csv(snakemake.output.solar)
    
    # wind data
    df = retrieve_wind_timeseries(regions, use_actual=use_actual)
    df = process_wind_timeseries(df)
    df.to_csv(snakemake.output.wind)