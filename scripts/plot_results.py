import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

    fig, ax = plt.subplots(figsize=(12, 6))

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
    ax.set_ylim(-n.storage_units_t.p_store.max().max()/1e3-2.5, n.loads_t.p_set.sum(axis=1).max()/1e3+2.5)
    plt.tight_layout()
    
    return fig, ax


def plot_capacity(n):
    
    fig, ax = plt.subplots(figsize=(12,8))
    
    df = pd.concat([n.generators[['p_nom', 'p_nom_opt']], 
                    n.storage_units[['p_nom','p_nom_opt']]])\
                        .replace(0, np.nan)\
                            .dropna(axis=0, how='all')
                            
    df.plot.bar(ax=ax)
    plt.tight_layout()
    
    return fig,ax
    
    
def plot_emissions(n):
    
    emissions_data = n.carriers\
                        .reset_index()\
                            .sort_values(by='Carrier')\
                                .set_index('Carrier')[['co2_emissions']]\
                                    .drop('Batteries')
    
    p_by_carrier = n.generators_t.p.groupby(n.generators.carrier, axis=1).sum()
    
    emissions_by_carrier_year = p_by_carrier.groupby(n.generators_t.p.index.get_level_values('timestep').year).sum()
    
    annual_emissions = emissions_by_carrier_year * emissions_data['co2_emissions'].sum(axis=1).to_frame()

    fig, ax = plt.subplots(figsize=(12,8))
    
    (annual_emissions/1e6).plot.bar(ax=ax)
    
    ax.set_ylabel("Mtons CO2/year")


if __name__ == "__main__":
    
    n = pypsa.Network(snakemake.input.solved_network)
    
    fig, ax = plot_dispatch(n, year=int(snakemake.config['plot_year']))
    
    plt.savefig(snakemake.output.dispatch_figure)
    
    fig, ax = plot_capacity(n)
    plt.savefig(snakemake.output.capacity_figure)