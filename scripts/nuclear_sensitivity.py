import pypsa
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
from glob import glob
from tqdm import tqdm
import seaborn as sb
from pathlib import Path

if __name__ == "__main__":
    results_path = Path("../pypsa-illinois/results/advanced_nuclear_sensitivity/")
    results_path.mkdir(parents=True, exist_ok=True)
    
    n = pypsa.Network("../pypsa-illinois/results/advanced_nuclear_v1.1/networks/illinois_solved.nc")
    
    delta = 0.25
    costs = np.arange(1, 3+delta, delta)
    
    smr_cost_2030 = n.generators.loc[n.generators.index.str.contains('SMR'), 'capital_cost'].unique()[0]
    
    results = []
    for cost in tqdm(costs):
        n_test = n.copy()
        n_test.generators.loc[n_test.generators.index.str.contains('SMR'), 'capital_cost'] = smr_cost_2030*cost
        
        n_test.optimize(solver_name='cplex')
        
        scenario = f"cost-2023_growth_0.01_demand-1.36E+08_atb-Moderate-X-{cost}_v1.1"
        
        results.append(n_test)
        
        network_folder = results_path/scenario
        network_folder.mkdir(exist_ok=True)
        
        n.export_to_netcdf(str(network_folder/"illinois_solved.nc"))