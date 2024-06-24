import networkx as nx
import pandas as pd
import numpy as np
import pypsa


def add_snapshots(network):
    resolution = int(snakemake.config['time_res'])

    model_years = np.array(snakemake.config['model_years']).astype('int')
    
    n_hours = 8760
    
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
    n.export_to_netcdf(snakemake.output.network)