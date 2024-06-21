import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from gridstatus import EIA
import os
from dotenv import load_dotenv
from pathlib import Path

curr_dir_os = Path(os.path.dirname(os.path.abspath(__file__)))

# load_dotenv("../.env")

if __name__ == "__main__":
    
    eia = EIA()
    route = "electricity/rto/region-sub-ba-data"
    
    start_date = '2019-01-01'
    end_date = '2024-01-01'
    
    facets = {
                'subba':snakemake.config['rto_subba'],  
            }
    
    demand = eia.get_dataset(dataset=route,
                start=start_date,
                end=end_date,
                facets=facets,
                n_workers=4,
                verbose=True)

    try:
        demand['Interval End'] = demand['Interval End'].dt.tz_localize(None)
    except:
        pass
    
    demand.MW = demand.MW.astype("float")
    
    demand_pivot = demand.pivot_table(index=['Interval End'],
                   columns='Subregion',
                   values='MW')
    
    demand_pivot[demand_pivot>60e3] = np.nan
    demand_pivot = demand_pivot.interpolate("linear")
    
    demand_pivot.rename(columns=dict(zip(snakemake.config['rto_subba'],
                                         snakemake.config['region_names'])),
                        inplace=True)
    
    demand_pivot.to_csv(snakemake.output.load)