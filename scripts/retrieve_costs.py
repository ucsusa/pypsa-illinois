import pandas as pd
from nrelpy.atb import ATBe


if __name__ == "__main__":
    atb_params = snakemake.config['atb_params']
    
    atb = ATBe(atb_params['atb_year'])
    df = atb.raw_dataframe

    df['technology_alias'] = df['technology_alias'].replace({"Utility-Scale Battery Storage":"Batteries",
                                                             "Biopower":"Biomass",
                                                             "Land-Based Wind":"Wind",
                                                             "Utility PV":"Solar"})


    df['techdetail'] = df['techdetail'].replace({'Nuclear':'LWR',
                                                 'Class1':'Land-Based Wind',
                                                 'Class2':'Utility PV',
                                                 'Dedicated':'Biopower'})

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
                                
    df_pivot.rename(columns={'Fixed O&M':'FOM',
                             'Variable O&M':'VOM',
                             },
                    inplace=True)
    # add natural gas and coal costs                            

    prices_url = "https://www.eia.gov/electricity/annual/html/epa_07_04.html"
    fuel_prices = pd.read_html(prices_url)[1].set_index(pd.Series(range(2012,2023))).T

    heatrate_url = "https://www.eia.gov/electricity/annual/html/epa_08_01.html"
    heatrates = pd.read_html(heatrate_url)[1].set_index("Year")

    price_col = 'Average Cost (Dollars per MMBtu)'

    coal_price = fuel_prices.loc[('Coal',
                                  'Subbituminous',
                                  price_col),
                                 atb_params['atb_year']]
                                
    naturalgas_price = fuel_prices.loc[('Natural Gas',
                                        slice(None), 
                                        price_col),
                                        atb_params['atb_year']].values[0]
    petroleum_price = fuel_prices.loc[('Petroleum',
                                       slice(None),
                                       price_col),
                                       atb_params['atb_year']].values[0]
    
    ng_heatrate = heatrates.at[atb_params['atb_year'], 'Natural Gas']/1e3  # converts kwh/btu to mwh/mmbtu
    coal_heatrate = heatrates.at[atb_params['atb_year'], 'Coal']/1e3
    petroleum_heatrate = heatrates.at[atb_params['atb_year'], 'Petroleum']/1e3
    
    df_pivot.loc[('Natural Gas',
                  slice(None),
                  slice(None)), 
                 'Fuel'] = naturalgas_price*ng_heatrate
    df_pivot.loc[('Coal',
                  slice(None),
                  slice(None)), 
                 'Fuel'] = coal_price*coal_heatrate
    # add petroleum
    df_t = df_pivot.T
    
    # these costs are from 2021 and should be updated to reflect inflation.
    df_t['Petroleum','Petroleum',snakemake.config['model_year']] = [1158, 
                                   27.94, 
                                   1.78, 
                                   petroleum_price*petroleum_heatrate]   # from here https://www.eia.gov/electricity/generatorcosts/
    df_pivot = df_t.T
    
    df_pivot.fillna(0., inplace=True)
    
    df_pivot.to_csv(snakemake.output.costs)