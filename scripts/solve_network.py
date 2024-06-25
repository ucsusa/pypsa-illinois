import pandas as pd
import numpy as np
import pypsa
import matplotlib.pyplot as plt


if __name__ == "__main__":
    
    n = pypsa.Network(snakemake.input.elec_network)
    
    n.optimize(solver_name='cplex')
    
    n.export_to_netcdf(snakemake.output.solved_network)