import pandas as pd
import numpy as np
import pypsa
import matplotlib.pyplot as plt

def plot_dispatch(n, year=2025, month=7):
    
    time = (year, f'{year}-0{month}')
    p_by_carrier = n.generators_t.p.groupby(n.generators.carrier, axis=1).sum().div(1e3)
    
    p_by_carrier = p_by_carrier[['Nuclear', 
                                 'Coal', 
                                 'Natural Gas', 
                                 'Biomass', 
                                 'Petroleum', 
                                 'Solar',
                                 'Wind']]

    if not n.storage_units.empty:
        sto = n.storage_units_t.p.T.groupby(n.storage_units.carrier).sum().T.div(1e3)
        p_by_carrier = pd.concat([p_by_carrier, sto], axis=1)

    fig, ax = plt.subplots(figsize=(6, 3))

    color = p_by_carrier.columns.map(n.carriers.color)

    p_by_carrier.where(p_by_carrier > 0).loc[time].plot.area(
        ax=ax,
        linewidth=0,
        color=color,
    )

    charge = p_by_carrier.where(p_by_carrier < 0).dropna(how="all", axis=1).loc[time]

    if not charge.empty:
        charge.plot.area(
            ax=ax,
            linewidth=0,
            color=charge.columns.map(n.carriers.color),
        )

    n.loads_t.p_set.sum(axis=1).loc[time].div(1e3).plot(ax=ax, c="k")

    ax.legend(loc=(1.05, 0))
    ax.set_ylabel("GW")
    ax.set_ylim(n.generators_t.p.min().min(), n.loads_t.p_set.sum(axis=1).max()/1e3+2.5)
    plt.tight_layout()
    
    return fig, ax


if __name__ == "__main__":
    
    n = pypsa.Network(snakemake.input.elec_network)
    
    n.optimize(solver_name='highs')
    
    n.export_to_netcdf(snakemake.output.solved_network)
    
    fig, ax = plot_dispatch(n, year=int(snakemake.config['plot_year']))
    
    plt.savefig(snakemake.output.dispatch_figure)