import networkx as nx
import pandas as pd
import numpy as np
import pypsa


model_years = np.array(snakemake.config['model_years']).astype('int')
resolution = int(snakemake.config['time_res'])
myopic = snakemake.config['myopic']

def add_snapshots(network, myopic=False):

    n_hours = 8760
    
    if myopic:
        raise NotImplementedError()
    else:
        snapshots = pd.DatetimeIndex([])
        for year in model_years:
            period = pd.date_range(start=f"{year}-01-01",
                                            freq=f"{resolution}h",
                                            periods=n_hours / resolution)
            snapshots = snapshots.append(period)
        
        network.snapshots = pd.MultiIndex.from_arrays([snapshots.year, snapshots])
        network.investment_periods = model_years
        
        network.snapshot_weightings.loc[:,:] = resolution
    
    return

def base_network():
    n = pypsa.Network()
    n.name = 'PyPSA-Illinois'
    
    buses = pd.read_csv(snakemake.input.buses, index_col='name')
    lines = pd.read_csv(snakemake.input.lines, index_col=0)
    
    line_config = snakemake.config['lines']
    v_nom = line_config['v_nom']

    
    n.import_components_from_dataframe(buses, 'Bus')
    n.import_components_from_dataframe(lines, 'Line')
    
    return n

if __name__ == "__main__":
    n = base_network()
    add_snapshots(n)
    
    year_diff = np.diff(model_years)
    if len(year_diff) == 0:
        pad = 1
    else:
        pad = year_diff[-1]
    
    n.investment_period_weightings["years"] = list(year_diff) + [pad]

    r = 0.01
    T = 0
    for period, nyears in n.investment_period_weightings.years.items():
        discounts = [(1 / (1 + r) ** t) for t in range(T, T + nyears)]
        n.investment_period_weightings.at[period, "objective"] = sum(discounts)
        T += nyears
    
    n.export_to_netcdf(snakemake.output.network)