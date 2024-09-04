import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

TECH_ORDER = ['Nuclear',
              'Coal',
              'Natural Gas',
              'Biomass',
              'Petroleum',
              'Solar',
              'Wind']


def power_by_carrier(n):
    p_by_carrier = n.generators_t.p.T.groupby(
        n.generators.carrier).sum().T 

    return p_by_carrier


def plot_dispatch(n, year=2025, month=7):

    time = (year, f'{year}-0{month}')
    p_by_carrier = power_by_carrier(n).div(1e3)

    p_by_carrier = p_by_carrier[TECH_ORDER]

    if not n.storage_units.empty:
        sto = n.storage_units_t.p.T.groupby(
            n.storage_units.carrier).sum().T.div(1e3)
        p_by_carrier = pd.concat([p_by_carrier, sto], axis=1)

    # y-limits
    y_min = -n.storage_units_t.p_store.max().max() / 1e3
    y_max = n.loads_t.p_set.sum(axis=1).max() / 1e3
    margin = 0.1
    y_low = (1 + margin) * y_min
    y_high = (1 + margin) * y_max

    fig, ax = plt.subplots(figsize=(12, 6))

    color = p_by_carrier.columns.map(n.carriers.color)

    p_by_carrier.where(p_by_carrier > 0).loc[time].plot.area(
        ax=ax,
        linewidth=0,
        color=color,
        ylim=(y_low - margin, y_high + margin)
    )

    charge = p_by_carrier.where(
        p_by_carrier < 0).dropna(
        how="all",
        axis=1).loc[time]

    if not charge.empty:
        charge.plot.area(
            ax=ax,
            linewidth=0,
            color=charge.columns.map(n.carriers.color),
            ylim=(y_low - margin, y_high + margin)
        )

    n.loads_t.p_set.sum(axis=1).loc[time].div(1e3).plot(ax=ax, c="k")

    ax.legend(loc=(1.05, 0))
    ax.set_ylabel("GW", fontsize=16)
    plt.tight_layout()

    return fig, ax


def plot_capacity(n):

    fig, ax = plt.subplots(figsize=(12, 8))

    df = pd.concat([n.generators[['p_nom', 'p_nom_opt']],
                    n.storage_units[['p_nom', 'p_nom_opt']]])\
        .replace(0, np.nan)\
        .dropna(axis=0, how='all')

    df.plot.bar(ax=ax)
    plt.tight_layout()

    return fig, ax


def plot_emissions(n, time_res):

    emissions_data = n.carriers\
        .reset_index()\
        .sort_values(by='Carrier')\
        .set_index('Carrier')[['co2_emissions']]\
        .drop('Batteries')

    emissions = n.generators_t.p * \
        n.generators.carrier.map(n.carriers.co2_emissions)
    total_emissions = n.snapshot_weightings.generators @ emissions.sum(
        axis=1).div(1e6)

    p_by_carrier = power_by_carrier(n)

    p_by_carrier_year = (
        p_by_carrier.groupby(
            p_by_carrier.index.get_level_values('period')).sum())

    annual_emissions = (
        p_by_carrier_year *
        emissions_data['co2_emissions']).sum(
        axis=1).to_frame()

    annual_emissions = annual_emissions * time_res / 1e6

    annual_emissions.columns = ['Mtonnes CO2/year']

    fig, ax = plt.subplots(figsize=(12, 8))

    annual_emissions.plot.bar(ax=ax)

    ax.set_ylabel("Mtonnes CO2/year")

    return fig, ax


def plot_active_units(n):

    fig, ax = plt.subplots(figsize=(12, 8))
    c = "StorageUnit"
    df = pd.concat(
        {
            period: n.get_active_assets(c, period) * n.df(c).p_nom_opt
            for period in n.investment_periods
        },
        axis=1,
    )
    df = df.groupby(n.storage_units.carrier).sum()

    c = "Generator"
    df2 = pd.concat(
        {
            period: n.get_active_assets(c, period) * n.df(c).p_nom_opt
            for period in n.investment_periods
        },
        axis=1,
    )

    df2 = df2.groupby(n.generators.carrier).sum()
    df = pd.concat([df, df2])
    color = df2.columns.map(n.carriers.color)
    df.T.plot.bar(
        ax=ax,
        stacked=True,
        edgecolor="white",
        width=1,
        ylabel="Capacity",
        xlabel="Investment Period",
        rot=0,
        figsize=(10, 5),
        # color=color
    )
    plt.tight_layout()

    return fig, ax


def plot_monthly_generation(n, time_res, model_years):
    p_by_carrier = power_by_carrier(n) * time_res

    p_by_carrier = p_by_carrier.resample('ME', level='timestep').sum()

    p_by_carrier = p_by_carrier[TECH_ORDER]

    color = p_by_carrier.columns.map(n.carriers.color)

    fig, ax = plt.subplots(1,len(model_years), figsize=(16, 8), sharey=True)
    for i, year in enumerate(model_years):
        p_by_carrier.loc[str(year)].plot.area(ax=ax[i],
                                color=color,
                                fontsize=16, )

        ax[i].set_xlabel('')

    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    return fig, ax


if __name__ == "__main__":

    n = pypsa.Network(snakemake.input.solved_network)

    fig, ax = plot_dispatch(n, year=int(snakemake.config['plot_year']))

    plt.savefig(snakemake.output.dispatch_figure)

    fig, ax = plot_capacity(n)
    plt.savefig(snakemake.output.capacity_figure)

    time_res = float(snakemake.config['time_res'])
    fig, ax = plot_emissions(n, time_res)
    plt.savefig(snakemake.output.emissions_figure)

    fig, ax = plot_active_units(n)
    plt.savefig(snakemake.output.active_units_figure)


    model_years = snakemake.config['model_years']
    fig, ax = plot_monthly_generation(n, time_res, model_years)
    plt.savefig(snakemake.output.monthly_generation_figure)
