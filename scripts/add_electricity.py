import networkx as nx
import pandas as pd
import numpy as np
import scipy as sp
import geopandas as gpd
from tqdm import tqdm
import pypsa

version = 22
model_years = np.array(snakemake.config['model_years']).astype('int')
resolution = int(snakemake.config['time_res'])
scale = snakemake.config['geo_res']
idx_opts = {"rto": "balancing_authority_code",
            "county": "county"}
growth_rates = snakemake.config['growth_rates']
pudl_year = int(snakemake.config['fuel_cost_year'])
wind_cf = float(snakemake.config['turbine_params']['capacity_factor'])
retirements_df = pd.DataFrame(snakemake.config['retirements']).fillna(0)
capacity_limits_df = pd.DataFrame(snakemake.config['capacity_max']).fillna(0)

BUILD_YEAR = 2025  # a universal build year place holder


def annuity(r, n):
    return r / (1 - 1 / (1 + r)**n)


def create_snapshots():
    timestamps = pd.DatetimeIndex([])
    for year in model_years:
        period = pd.date_range(start=f"{year}-01-01",
                               freq=f"{resolution}h",
                               periods=8760 / resolution)
        timestamps = timestamps.append(period)

    return timestamps


def load_costs():

    costs = pd.read_csv(
        snakemake.input.costs,
        index_col=[
            'technology_alias',
            'techdetail',
            'core_metric_variable'])

    costs[['CAPEX', 'FOM']] *= 1e3  # convert from /kW to /MW

    r = float(snakemake.config['discount_rate'])
    lifetimes = snakemake.config['lifetime']

    costs = costs.assign(rate=r, lifetime=20)

    # assign lifetimes to technologies
    carriers = snakemake.config['atb_params']['carrier']
    for carrier in carriers:
        costs.loc[(carrier, slice(None), slice(None)),
                  'lifetime'] = float(lifetimes[carrier])

    annuity_col = annuity(costs["rate"], costs["lifetime"])

    costs = costs.assign(
        capital_cost=(
            annuity_col *
            costs['CAPEX'] +
            costs['FOM']),
        marginal_cost=(
            costs['VOM'] +
            costs['Fuel']))

    return costs


def load_costs_ts():

    fuel_costs = pd.read_csv(snakemake.input.fuel_cost_timeseries,
                             parse_dates=True,
                             index_col=['report_date'])

    return fuel_costs


def load_existing_generators():
    generators = pd.read_csv(snakemake.input.generators,
                             index_col=idx_opts[scale])

    return generators


def load_build_years():
    build_years = pd.read_csv(snakemake.input.build_years,
                              index_col=idx_opts[scale])

    return build_years


def linear_growth(init_value, start_year, growth_rate, end_year=2050):
    def model(x, init_val, start, rate):
        return rate * init_val * (x - start) + init_val
    years = np.arange(start_year, end_year, 1).astype('int')
    growth_data = model(years, init_value, start_year, growth_rate)

    growth_df = pd.DataFrame({'demand': growth_data})
    growth_df.index = pd.date_range(start=str(start_year),
                                    periods=(end_year - start_year),
                                    freq='YE')

    return growth_df


# attach components
def attach_load(n):
    initial_demand = float(snakemake.config['total_demand'])

    if len(model_years) == 1:
        demand = [initial_demand]
    else:
        rate = float(snakemake.config['load_growth'])
        start_year = model_years[0]
        end_year = model_years[-1]
        N_years = end_year - start_year
        demand = linear_growth(initial_demand,
                               model_years[0],
                               rate
                               )
        demand = demand.loc[demand.index.year.isin(
            model_years)].values.flatten()

    load = pd.read_csv(
        snakemake.input.load,
        parse_dates=True,
        index_col="Interval End")

    start = int(snakemake.config['load_year'])
    n_model_years = len(snakemake.config['model_years'])
    end = start + n_model_years - 1
    load = load.loc[str(start):str(end)][:int(n_model_years * 8760)]
    # normalize the load data

    for d, y in zip(demand, load.index.year.unique()):
        load.loc[str(y)] = load.loc[str(y)].div(
            load.loc[str(y)].sum(axis=0).sum(), axis=1) * d

    load = load.resample(f"{snakemake.config['time_res']}h").mean()
    load = load[:len(n.snapshots)]
    load.set_index(n.snapshots, inplace=True)

    print('Adding loads to model')
    for bus in tqdm(n.buses.index):
        n.add(class_name="Load",
              name=bus,
              bus=bus,
              p_set=load[bus]
              )

    return


def load_emissions():
    df = pd.read_csv(snakemake.input.emissions, index_col=0)

    return df


def add_carriers(n, costs, emissions):
    available_carriers = list(costs['technology_alias'].unique())

    for carrier in available_carriers:
        # emissions data
        try:
            co2_emissions = emissions.at[carrier, 'tCO2 per MWh']
        except KeyError:
            co2_emissions = 0.00

        n.add(class_name="Carrier",
              name=carrier,
              co2_emissions=co2_emissions,
              color=snakemake.config['carrier_colors'][carrier],
              max_growth=float(growth_rates[carrier]))
    return


def load_re_profile(n, carrier='Solar'):
    file_name = {'Solar': snakemake.input.solar_profile,
                 'Wind': snakemake.input.wind_profile}
    re_profile = pd.read_csv(file_name[carrier], parse_dates=True, index_col=0)
    re_profile = re_profile.resample(
        f"{snakemake.config['time_res']}h").mean().dropna(axis=0)
    re_profile.set_index(n.snapshots, inplace=True)

    return re_profile


def attach_renewables(
        n,
        costs,
        generators=None,
        build_years=None,
        model_year=None):
    carriers = ['Wind', 'Solar']

    for bus in n.buses.index:
        for item in costs.itertuples():
            tech = item.techdetail
            carrier = item.technology_alias
            if carrier in carriers:
                re_profile = load_re_profile(n, carrier=carrier)
                # existing capacity
                if isinstance(
                        generators,
                        pd.DataFrame) and isinstance(
                        build_years,
                        pd.DataFrame):
                    try:
                        p_nom = generators.at[bus, tech]
                        build_year = build_years.at[bus, tech]
                    except BaseException:
                        continue
                    name = f"{bus} {tech} EXIST"
                    extendable = False
                elif model_year:
                    build_year = model_year
                    p_nom = 0.0
                    name = f"{bus} {tech} {model_year}"
                    extendable = tech in snakemake.config['extendable_techs']
                    build_year = model_year

                p_max_pu = re_profile[bus]

                n.add(class_name="Generator",
                      name=name,
                      bus=bus,
                      p_nom=p_nom,
                      p_nom_min=p_nom,
                      p_max_pu=p_max_pu,
                      p_nom_extendable=extendable,
                      carrier=carrier,
                      capital_cost=item.capital_cost,
                      marginal_cost=item.marginal_cost,
                      lifetime=item.lifetime,
                      build_year=build_year)
    return


def attach_generators(
        n,
        costs,
        generators=None,
        build_years=None,
        model_year=None,
        costs_ts=None):
    carriers = ['Natural Gas', 'Biomass', 'Coal', 'Nuclear', 'Petroleum']
    for bus in n.buses.index:
        for item in costs.itertuples():
            tech = item.techdetail
            carrier = item.technology_alias
            if carrier in carriers:
                # existing capacity
                if isinstance(
                        generators,
                        pd.DataFrame) and isinstance(
                        build_years,
                        pd.DataFrame):
                    try:
                        p_nom = generators.at[bus, tech]
                        build_year = build_years.at[bus, tech]
                    except BaseException:
                        continue
                    name = f"{bus} {tech} EXIST"
                    extendable = False
                elif model_year:
                    build_year = model_year
                    p_nom = 0.0
                    name = f"{bus} {tech} {model_year}"
                    extendable = tech in snakemake.config['extendable_techs']
                    build_year = model_year

                # ramp limits
                if tech in ['LWR']:
                    ramp_limit_up = 0.01 * resolution
                    ramp_limit_down = 0.01 * resolution
                    p_nom_min = p_nom
                elif tech in ['IGCCAvgCF', 'Biopower']:
                    ramp_limit_up = 0.1 * resolution
                    ramp_limit_down = 0.1 * resolution
                    p_nom_min = 0
                elif tech in ['CCAvgCF', 'CTAvgCF']:
                    ramp_limit_up = min(0.6 * resolution, 1.0)
                    ramp_limit_down = min(0.6 * resolution, 1.0)
                    p_nom_min = 0
                elif tech in ['Petroleum']:
                    ramp_limit_up = 1.0
                    ramp_limit_down = 1.0
                    p_nom_min = 0
                else:
                    ramp_limit_up = 1.0
                    ramp_limit_down = 1.0
                    p_nom_min = p_nom

                # minimum/maximum power output
                if tech == 'LWR':
                    p_min_pu = 0.0
                    p_max_pu = 1.0
                else:
                    p_max_pu = 1
                    p_min_pu = 0

                # time series marginal costs
                if ((tech in ['CTAvgCF', 'CCAvgCF', 'IGCCAvgCF'])
                        and isinstance(costs_ts, pd.DataFrame)):

                    # select year to replicate
                    cost_data = costs_ts.loc[str(pudl_year), carrier].values

                    # get the VOM cost
                    cost_data = cost_data + item.VOM

                    marginal_cost = np.tile(cost_data, len(model_years))
                else:
                    marginal_cost = item.marginal_cost

                n.add(class_name="Generator",
                      name=name,
                      bus=bus,
                      p_nom=p_nom,
                      p_nom_min=p_nom_min,
                      p_nom_extendable=extendable,
                      p_min_pu=p_min_pu,
                      p_max_pu=p_max_pu,
                      carrier=carrier,
                      capital_cost=item.capital_cost,
                      marginal_cost=marginal_cost,
                      lifetime=item.lifetime,
                      ramp_limit_down=ramp_limit_down,
                      ramp_limit_up=ramp_limit_up,
                      build_year=build_year,
                      )
    return


def attach_storage(
        n,
        costs,
        generators=None,
        build_years=None,
        model_year=None):
    carriers = ["Batteries"]
    for bus in n.buses.index:
        for item in costs.itertuples():
            tech = item.techdetail
            carrier = item.technology_alias
            if carrier in carriers:
                # existing capacity
                if isinstance(
                        generators,
                        pd.DataFrame) and isinstance(
                        build_years,
                        pd.DataFrame):
                    try:
                        p_nom = generators.at[bus, tech]
                        build_year = build_years.at[bus, tech]
                    except BaseException:
                        continue
                    name = f"{bus} {tech} EXIST"
                    extendable = False
                elif model_year:
                    build_year = model_year
                    p_nom = 0.0
                    name = f"{bus} {tech} {model_year}"
                    extendable = tech in snakemake.config['extendable_techs']
                    build_year = model_year

                n.add(class_name="StorageUnit",
                      name=name,
                      bus=bus,
                      p_nom=p_nom,
                      p_nom_min=p_nom,
                      p_nom_extendable=extendable,
                      carrier=carrier,
                      capital_cost=item.capital_cost,
                      marginal_cost=item.marginal_cost,
                      lifetime=item.lifetime,
                      max_hours=float(tech.split(' ')[0].strip('Hr')),
                      cyclic_state_of_charge=False,
                      build_year=build_year)

    return


def add_retirements(n):
    
    carriers = retirements_df.columns
    
    df = retirements_df.loc[model_years, :].cumsum()
    c_by_carrier = n.generators.groupby('carrier').sum().loc[carriers, 'p_nom']
    
    for carrier in carriers:
        existing_cap = c_by_carrier.loc[carrier]
        for year in model_years:
            try:
                retirement = df.loc[year, carrier]
                remaining = max(existing_cap - retirement, 0)
                
                limit = remaining * (8760/resolution)
                
                n.add(class_name="GlobalConstraint",
                        name=f'{carrier} Energy Limit {year}',
                        sense='<=',
                        carrier_attribute=carrier,
                        constant=limit,
                        type='operational_limit',
                        investment_period=year,
                        )
            except (AttributeError, KeyError, TypeError):
                pass
            
    return


def add_capacity_max(n):
    
    carriers = capacity_limits_df.columns
    
    df = capacity_limits_df.copy()
    n_buses = n.buses.shape[0]
    
    for carrier in carriers:
        for year in df.index:
            try:
                limit = df.loc[year, carrier] / n_buses
                n.add(class_name="GlobalConstraint",
                        name=f'{carrier} Capacity Limit {year}',
                        sense='<=',
                        carrier_attribute=carrier,
                        constant=limit,
                        type='tech_capacity_expansion_limit',
                        investment_period=year,
                        )
            except (AttributeError, KeyError, TypeError):
                pass
            
    return
    

if __name__ == "__main__":

    n = pypsa.Network(snakemake.input.base_network)
    costs = load_costs()
    costs.to_csv("data/final_costs.csv")
    current_costs = costs.xs((slice(None), slice(None), 2020))
    current_costs = current_costs.reset_index()
    current_costs.loc[current_costs['technology_alias']
                      == 'Solar', 'techdetail'] = 'Utility PV'

    generators = load_existing_generators()
    build_years = load_build_years()
    emissions = load_emissions()
    costs_ts = load_costs_ts()

    attach_load(n)
    add_carriers(n,
                 costs=current_costs,
                 emissions=emissions)

    # add existing technology
    attach_renewables(n,
                      costs=current_costs,
                      generators=generators,
                      build_years=build_years
                      )
    attach_generators(n,
                      costs=current_costs,
                      generators=generators,
                      build_years=build_years,
                      costs_ts=costs_ts
                      )
    attach_storage(n,
                   costs=current_costs,
                   generators=generators,
                   build_years=build_years
                   )

    # add new technology
    for year in model_years:
        # current_costs = costs.xs((slice(None), slice(None), year))
        # current_costs = current_costs.reset_index()
        # current_costs.loc[current_costs['technology_alias']=='Solar', 'techdetail'] = 'Utility PV'
        attach_renewables(n,
                          costs=current_costs,
                          model_year=year
                          )
        attach_generators(n,
                          costs=current_costs,
                          model_year=year,
                          costs_ts=costs_ts
                          )
        attach_storage(n,
                       costs=current_costs,
                       model_year=year
                       )


    # add energy constraints
    add_retirements(n)
    
    # add capacity constraints
    add_capacity_max(n)

    # add co2 constraint
    emissions_dict = snakemake.config['co2_limits']
    for y in model_years:
        # if y in emissions_dict.keys():
        try:
            n.add(class_name="GlobalConstraint",
                  name=f"CO2 Limit {y}",
                  carrier_attribute='co2_emissions',
                  sense="<=",
                  investment_period=y,
                  constant=float(emissions_dict[y]) * 1e6)
        # else:
        except (AttributeError, KeyError, TypeError):
            pass

    # modify wind capacity factor
    wind_gen = n.generators[n.generators.carrier == 'Wind'].index
    n.generators_t.p_max_pu.loc[:, wind_gen] = ((n.generators_t.p_max_pu[wind_gen] / (
        n.generators_t.p_max_pu[wind_gen].sum() / (len(n.snapshots))) * wind_cf))

    n.export_to_netcdf(snakemake.output.elec_network)
