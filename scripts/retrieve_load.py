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
                'parent':['MISO','PJM']
                # 'subba':['0004','CE'],  
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
    
    demand_pivot = demand.pivot_table(index=['Interval End','BA'],
                   columns='Subregion',
                   values='MW')
    
    
    pjm_demand = demand_pivot.xs((slice(None), 'PJM')).dropna(axis=1)
    miso_demand = demand_pivot.xs((slice(None), 'MISO')).dropna(axis=1)
    
    pjm_demand = pjm_demand.reset_index().set_index("Interval End")
    miso_demand = miso_demand.reset_index().set_index("Interval End")
    
    pjm_demand.to_csv(str(curr_dir_os/"../data/pjm_demand.csv"))
    miso_demand.to_csv(str(curr_dir_os/"../data/miso_demand.csv"))