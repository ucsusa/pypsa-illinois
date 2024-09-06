import pandas as pd
import numpy as np
import pypsa
import matplotlib.pyplot as plt


if __name__ == "__main__":
    solver = snakemake.config['solver']
    multi_invest = snakemake.config['multi_investment_periods']

    n = pypsa.Network(snakemake.input.elec_network)

    try:
        n.optimize(solver_name=solver, multi_investment_periods=multi_invest)
    except:
        print('No CPLEX solver found. Reverting to HiGHS')
        n.optimize(solver_name='highs', multi_investment_periods=multi_invest)

    n.export_to_netcdf(snakemake.output.solved_network)