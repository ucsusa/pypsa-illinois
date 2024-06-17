import pandas as pd
from nrelpy.atb import ATBe


if __name__ == "__main__":
    atb_params = snakemake.config['atb_params']
    
    atb = ATBe(atb_params['atb_year'])
    df = atb.raw_dataframe

    df_pivot = df.pivot_table(index=['scenario',
                                    'core_metric_case',
                                    'crpyears',
                                    'technology_alias',
                                    'techdetail',
                                    'core_metric_variable'],
                              columns='core_metric_parameter',
                              values='value')
    df_pivot = df_pivot.xs((atb_params['scenario'],
                            atb_params['case'],
                            str(atb_params['crp']),
                            slice(None),
                            slice(None),
                            slice(None)
                            ))
    
    df_pivot = df_pivot.loc[(atb_params['carrier'],
                             atb_params['technology'] ,
                             slice(None)), 
                            ['CAPEX', 'Fixed O&M',
                             'Variable O&M','Fuel']]\
                                .dropna(axis=0, how='all')\
                                .drop_duplicates(keep='first')
                                
    # add natural gas and coal costs                            

    url = "https://www.eia.gov/electricity/annual/html/epa_07_04.html"
    fuel_prices = pd.read_html(url)[1].set_index(pd.Series(range(2012,2023))).T

    coal_price = fuel_prices.loc[('Coal',
                                  'Bituminous',
                                  'Average Cost (Dollars per MMBtu)'),
                                 atb_params['atb_year']]
                                
    naturalgas_price = fuel_prices.loc[('Natural Gas',
                                  slice(None), 
                                  'Average Cost (Dollars per MMBtu)'),
                                 2022].values[0]
    
    ng_heatrate = snakemake.config['tech_params']['natural_gas']['heat_rate']
    coal_heatrate = snakemake.config['tech_params']['coal']['heat_rate']
    
    df_pivot.loc[('Natural Gas',
                  slice(None),
                  slice(None)), 
                 'Fuel'] = naturalgas_price*ng_heatrate
    df_pivot.loc[('Coal',
                  slice(None),
                  slice(None)), 
                 'Fuel'] = naturalgas_price*ng_heatrate
    
    df_pivot.to_csv(snakemake.output.costs)