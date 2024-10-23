import numpy as np
import pandas as pd
from nrelpy.atb import ATBe
import sys

sys.path.append("functions")
from data_functions import load_heatrates


n_illinois_reactors = 11
total_lwr_capacity = 12415.1
average_reactor_size = total_lwr_capacity / n_illinois_reactors  # MW
nuclear_params = snakemake.config['nuclear_params']
capacity_factor = float(nuclear_params['capacity_factor'])
capital_cost = float(nuclear_params['capital_cost'])  # $/MWh
nuclear_fuel = float(nuclear_params['fuel'])  # $/MWh
fixed_om = float(nuclear_params['fixed_om'])  # $/MWh

# convert capital and fom to $/MW/yr
operating_hours = 8760 * capacity_factor
nuclear_cap_cost = capital_cost * operating_hours / 1e3  # $/kW/yr
nuclear_fixed_om = fixed_om * operating_hours / 1e3  # $/kW/yr


if __name__ == "__main__":
    atb_params = snakemake.config['atb_params']

    atb = ATBe(atb_params['atb_year'])
    df = atb.raw_dataframe

    df['technology_alias'] = df['technology_alias'].replace(
        {
            "Utility-Scale Battery Storage": "Batteries",
            "Biopower": "Biomass",
            "Land-Based Wind": "Wind",
            "Utility PV": "Solar"})

    df['techdetail'] = df['techdetail'].replace({'Nuclear': 'LWR',
                                                 'Class1': 'Land-Based Wind',
                                                 'Class2': 'Utility PV',
                                                 'Dedicated': 'Biopower'})

    df_pivot = df.pivot_table(index=['scenario',
                                     'core_metric_case',
                                     'crpyears',
                                     'technology_alias',
                                     'techdetail',
                                     'core_metric_variable'],
                              columns='core_metric_parameter',
                              values='value')
    # breakpoint()
    df_pivot = df_pivot.xs((atb_params['scenario'], 
                            atb_params['case'], 
                            slice(None), 
                            slice(None), 
                            slice(None), 
                            slice(None)))
    
    wacc_data = df_pivot.loc[('*', 
                              atb_params['carrier'], 
                              slice(None), 
                              slice(None)), 
                             ['WACC Nominal', 'WACC Real']]
    wacc_data.to_csv("data/wacc.csv")
    
    # df_pivot = df_pivot.xs((atb_params['scenario'],
    #                         atb_params['case'],
    #                         str(atb_params['crp']),
    #                         slice(None),
    #                         slice(None),
    #                         slice(None)
    #                         ))
    
    df_pivot = df_pivot.loc[(str(atb_params['crp']),
                             atb_params['carrier'],
                             atb_params['technology'],
                             slice(None)),
                            ['CAPEX', 'Fixed O&M',
                             'Variable O&M', 'Fuel']]

    df_pivot.rename(columns={'Fixed O&M': 'FOM',
                             'Variable O&M': 'VOM',
                             },
                    inplace=True)
    
    df_pivot = df_pivot.reset_index()
    wacc_data = wacc_data.reset_index()
    
    wacc_data['crpyears'] = str(atb_params['crp'])
    
    merged_data = df_pivot.merge(wacc_data, 
                                 on=['technology_alias',
                                     'core_metric_variable',
                                     'crpyears'],
                                 how='left')
    
    merged_data = merged_data.drop(columns=['techdetail_y','crpyears'])
    merged_data.rename(columns={'techdetail_x':'techdetail'}, inplace=True)
    merged_data = merged_data.loc[~((merged_data['technology_alias']=='Wind') 
                                    & (merged_data['techdetail'] == 'Utility PV'))]
    merged_data = merged_data.loc[~((merged_data['technology_alias']=='Solar') 
                                    & (merged_data['techdetail'] == 'Land-Based Wind'))]
    
    df_pivot = merged_data.set_index(['technology_alias','techdetail','core_metric_variable'])
    
    
    wacc_data = wacc_data.drop(columns=['crpyears', 
                                        'techdetail']).set_index(['technology_alias',
                                                                  'core_metric_variable'])
    # add natural gas and coal costs

    prices_url = "https://www.eia.gov/electricity/annual/html/epa_07_04.html"
    fuel_prices = pd.read_html(prices_url)[1].set_index(
        pd.Series(range(2012, 2023))).T

    heatrates = load_heatrates()

    price_col = 'Average Cost (Dollars per MMBtu)'

    coal_price = fuel_prices.loc[('Coal',
                                  'Subbituminous',
                                  price_col),
                                 atb_params['atb_year']]

    naturalgas_price = fuel_prices.loc[('Natural Gas',
                                        slice(None),
                                        price_col),
                                       atb_params['atb_year']].values[0]
    # $/mmbtu, 2023 avg https://www.eia.gov/todayinenergy/detail.php?id=61183#
    naturalgas_price = 2.57
    # naturalgas_price =

    petroleum_price = fuel_prices.loc[('Petroleum',
                                       slice(None),
                                       price_col),
                                      atb_params['atb_year']].values[0]

    # converts btu/kwh to mmbtu/mwh
    ng_heatrate = heatrates.at[atb_params['atb_year'], 'Natural Gas'] / 1e3
    coal_heatrate = heatrates.at[atb_params['atb_year'], 'Coal'] / 1e3
    petroleum_heatrate = heatrates.at[atb_params['atb_year'],
                                      'Petroleum'] / 1e3

    df_pivot.loc[('Natural Gas',
                  slice(None),
                  slice(None)),
                 'Fuel'] = naturalgas_price * ng_heatrate
    df_pivot.loc[('Coal',
                  slice(None),
                  slice(None)),
                 'Fuel'] = coal_price * coal_heatrate
    # add petroleum
    df_t = df_pivot.T
    # [capital cost, fixed om cost, variable om cost, fuel cost]
    for year in range(2020, 2051, 1):
        # update costs for existing nuclear, these costs are from 2022
        df_t['Nuclear', 'LWR', year] = [
            nuclear_cap_cost, 
            nuclear_fixed_om, 
            0, 
            nuclear_fuel, 
            wacc_data.loc[('Nuclear',year), 'WACC Nominal'],
            wacc_data.loc[('Nuclear',year), 'WACC Real'],]

        # these costs are from 2021 and should be updated to reflect inflation.
        # from here https://www.eia.gov/electricity/generatorcosts/
        df_t['Petroleum', 'Petroleum', year] = [1158,
                                                27.94,
                                                1.78,
                                                (petroleum_price *
                                                 petroleum_heatrate),
                                                float(snakemake.config['discount_rate']),
                                                float(snakemake.config['discount_rate'])]
    df_pivot = df_t.T

    # set values for the battery WACC
    df_pivot.loc[('Batteries', 
                '4Hr Battery Storage', 
                slice(None)), 
                ['WACC Real']] = float(snakemake.config['battery_wacc_real'])
    df_pivot.loc[('Batteries', 
                '4Hr Battery Storage', 
                slice(None)), 
                ['WACC Nominal']] = float(snakemake.config['battery_wacc_nom'])

    df_pivot.fillna(0., inplace=True)

    df_pivot.to_csv(snakemake.output.costs)
    heatrates.to_csv(snakemake.output.heatrates)
