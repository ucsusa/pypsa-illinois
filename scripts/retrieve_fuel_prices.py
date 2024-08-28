import pandas as pd
import matplotlib.pyplot as plt
import requests
import numpy as np
import json
import os
from dotenv import load_dotenv
load_dotenv("../.env")


def retrieve_eia923m_pudl(state):
    url = (f"https://data.catalyst.coop/pudl/"
           f"out_eia923__monthly_fuel_receipts_costs.csv?"
           f"_labels=on&_stream=on&_sort=rowid&"
           f"state__exact={state}&_size=max")
    
    df = pd.read_csv(url)
    
    return df


def fuel_cost_pivot(dataframe):
    cost_pivot = dataframe.pivot_table(index=['report_date'],
                            columns=['fuel_type_code_pudl'],
                            values=['fuel_cost_per_mmbtu'])
    
    # removes first column index, leaves the fuel name
    cost_pivot = cost_pivot.droplevel(level=0, axis=1)
    
    cost_pivot = cost_pivot.rename({"coal":"Coal", 'gas':"Natural Gas", 'oil':'Petroleum'})
    
    return cost_pivot

def prepare_cost_timeseries(df, time_res):
    
    dataframe = df.copy()
    
    dataframe.index = pd.to_datetime(dataframe.index)
    # remove NaN values
    dataframe = dataframe.interpolate()
    
    dataframe = dataframe.resample(f'{time_res}h').ffill()
    
    return dataframe

if __name__ == "__main__":
    
    state_abbr = snakemake.config['state_abbr']
    time_res = snakemake.config['time_res']
    
    eia_923m_df = retrieve_eia923m_pudl(state=state_abbr)
    
    fuel_cost_df = fuel_cost_pivot(eia_923m_df)
    
    fuel_cost_df.to_csv(snakemake.output.fuel_costs)
    
    fuel_cost_resampled = prepare_cost_timeseries(fuel_cost_df, time_res)
    fuel_cost_resampled.to_csv(snakemake.output.fuel_cost_timeseries)