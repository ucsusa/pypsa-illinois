import networkx as nx
import pandas as pd
import numpy as np
import pypsa


def base_network():
    n = pypsa.Network()
    n.name = 'PyPSA-Illinois'
    
    buses = pd.read_csv(snakemake.input.buses, index_col='name')
    lines = pd.read_csv(snakemake.input.lines, index_col=0)
    
    line_config = snakemake.config['lines']
    v_nom = line_config['v_nom']

    model_year = snakemake.config['model_year']
    n.set_snapshots(pd.date_range(start=f"{model_year}-01-01", 
                                  end=f"{model_year+1}-01-01", 
                                  inclusive='left',
                                  freq=f"{snakemake.config['time_res']}h"))
    n.import_components_from_dataframe(buses, 'Bus')
    n.import_components_from_dataframe(lines, 'Line')
    
    return n

if __name__ == "__main__":
    n = base_network()
    n.export_to_netcdf(snakemake.output.network)