import numpy as np
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from us import states
import geopandas as gpd
from shapely import Point
from pathlib import Path
from tqdm import tqdm

sys.path.append("functions")

from nrel_data_api import parameters, make_csv_url

def handle_datetime(dataframe):
    """
    Combines time columns into a single timestamp column.
    Expects columns ['year','month','day','hour'].

    Parameters
    ----------
    dataframe : :class:`pd.DataFrame`
        A pandas dataframe.
    """
    frame = dataframe.copy()
    model_year = snakemake.config['model_year']
    timestamps = pd.date_range(f"{model_year}-01-01",f"{model_year+1}-01-01",
                               freq='1h', inclusive='left')
    
    try:
        frame.set_index(timestamps,inplace=True)
        frame.drop(columns=['Year','Month','Day','Hour','Minute'],inplace=True)
    except:
        raise ValueError
    
    return frame
    
    

def retrieve_solar_timeseries(region):
    """
    Retrieves data from NREL's national solar radiation database (NSRDB).

    Parameters
    ----------
    region : :class:`gpd.GeoDataFrame`
        A geopandas dataframe containing modeled bus regions.
    """
    parameters['attr_list'] = ['ghi']
    parameters['year'] = int(snakemake.config['solar_year'])
    pbar = tqdm(region[['name','x','y']].values, position=0, leave=True)
    frames = []
    for n, i, j in pbar:
        pbar.set_description(f"Processing {n}")
        parameters['lon'] = i
        parameters['lat'] = j
        URL = make_csv_url(parameters=parameters, 
                           kind='solar')
        df = pd.read_csv(URL, skiprows=2)
        df.rename(columns={'GHI':f"{n}"}, inplace=True)
        df = handle_datetime(df)
        frames.append(df)
    
    solar_df = pd.concat(frames, axis=1)
    
    return solar_df

def retrieve_wind_timeseries(region):
    """
    Retrieves data from NREL's national solar radiation database (NSRDB).

    Parameters
    ----------
    region : :class:`gpd.GeoDataFrame`
        A geopandas dataframe containing modeled bus regions.
    """
    wind_attr = ['windspeed_80m']
    parameters['attr_list'] = wind_attr
    parameters['year'] = int(snakemake.config['wind_year'])
    pbar = tqdm(region[['name','x','y']].values, position=0, leave=True)
    frames = []
    for n, i, j in pbar:
        pbar.set_description(f"Processing {n}")
        parameters['lon'] = i
        parameters['lat'] = j
        URL = make_csv_url(parameters=parameters, 
                           kind='wind')
        df = pd.read_csv(URL, skiprows=1)
        df.rename(columns={'wind speed at 80m (m/s)':f"{n}"}, inplace=True)
        df = handle_datetime(df)
        frames.append(df)
    
    wind_df = pd.concat(frames, axis=1)
    
    return wind_df


if __name__ == "__main__":
    
    regions = gpd.read_file(snakemake.input.supply_regions)
    
    df = retrieve_solar_timeseries(regions)
    
    df.to_csv(snakemake.output.solar)
    
    df = retrieve_wind_timeseries(regions)
    df.to_csv(snakemake.output.wind)