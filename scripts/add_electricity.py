import networkx as nx
import pandas as pd
import numpy as np
import scipy as sp
import geopandas as gpd
from tqdm import tqdm
import pypsa

model_years = np.array(snakemake.config['model_years']).astype('int')
resolution = int(snakemake.config['time_res'])


def annuity(r, n):
    return r / (1-1/(1+r)**n)


def create_snapshots():
    timestamps = pd.DatetimeIndex([])
    for year in model_years:
        period = pd.date_range(start=f"{year}-01-01",
                                        freq=f"{resolution}h",
                                        periods=8760/resolution)
        timestamps = timestamps.append(period)
    
    return timestamps


def load_costs():
    
    costs = pd.read_csv(snakemake.input.costs, index_col=['technology_alias',
                                                          'techdetail',
                                                          'core_metric_variable'])
    
    costs[['CAPEX','FOM']] *= 1e3  # convert from /kW to /MW
    
    r = float(snakemake.config['discount_rate'])
    lifetimes = snakemake.config['lifetime']
    
    costs = costs.assign(rate=r,lifetime=20)
    
    # assign lifetimes to technologies
    carriers = snakemake.config['atb_params']['carrier']
    for carrier in carriers:
        costs.loc[(carrier,slice(None), slice(None)), 'lifetime'] = float(lifetimes[carrier])
    
    annuity_col = annuity(costs["rate"], costs["lifetime"])
    
    costs = costs.assign(capital_cost = (annuity_col*costs['CAPEX'] + costs['FOM']),
                         marginal_cost = (costs['VOM']+costs['Fuel']))
    
    return costs

def load_existing_generators():
    generators = pd.read_csv(snakemake.input.generators, 
                             index_col='balancing_authority_code')

    return generators
    

# attach components
def attach_load(n):
    load = pd.read_csv(snakemake.input.load, parse_dates=True, index_col="Interval End")
    
    start = int(snakemake.config['load_year'])
    n_model_years = len(snakemake.config['model_years'])
    end = start + n_model_years - 1
    load = load.loc[str(start):str(end)][:int(n_model_years*8760)]
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
            co2_emissions = emissions.at[carrier,'tCO2 per MWh']
        except KeyError:
            co2_emissions = 0.0
            
        n.add(class_name="Carrier",
              name=carrier, 
              co2_emissions=co2_emissions,
              color=snakemake.config['carrier_colors'][carrier])
    return

def attach_renewables(n, costs, generators):
    carriers = ['Wind','Solar']
    wind_profile = pd.read_csv(snakemake.input.wind_profile, parse_dates=True, index_col=0)
    wind_profile = wind_profile.resample(f"{snakemake.config['time_res']}h").mean()
    wind_profile = pd.concat([wind_profile.loc[str(year)] for year in model_years])
    wind_profile = wind_profile[:len(n.snapshots)]
    wind_profile.set_index(n.snapshots, inplace=True)
    
    solar_profile = pd.read_csv(snakemake.input.solar_profile, parse_dates=True, index_col=0)
    solar_profile = solar_profile.resample(f"{snakemake.config['time_res']}h").mean()
    solar_profile = pd.concat([solar_profile.loc[str(year)] for year in model_years])
    solar_profile = solar_profile[:len(n.snapshots)]
    solar_profile.set_index(n.snapshots, inplace=True) 
     
    for bus in n.buses.index:
        for item in costs.itertuples():
            tech = item.techdetail
            carrier = item.technology_alias
            if carrier in carriers:
                # existing capacity
                if tech in generators.loc[bus].index:
                    p_nom = generators.at[bus, tech]
                else:
                    p_nom = 0.0
                    
                # renewables
                if carrier == 'Wind':
                    p_max_pu = wind_profile[bus]
                elif carrier == 'Solar':
                    p_max_pu = solar_profile[bus]
                    
                extendable = tech in snakemake.config['extendable_techs']

                n.add(class_name="Generator",
                    name=f"{bus} {tech}",
                    bus=bus,
                    p_nom=p_nom,
                    p_nom_min=p_nom,
                    p_max_pu=p_max_pu,
                    p_nom_extendable=extendable,
                    carrier=carrier,
                    capital_cost=item.capital_cost,
                    marginal_cost=item.marginal_cost,
                    lifetime=item.lifetime,)
    return 


def attach_generators(n, costs, generators):
    carriers = ['Natural Gas','Biomass','Coal','Nuclear','Petroleum']
    for bus in n.buses.index:
        for item in costs.itertuples():
            tech = item.techdetail
            carrier = item.technology_alias
            if carrier in carriers:
                # ramp limits
                if tech in ['LWR']:
                    ramp_limit_up = 0.01*resolution
                    ramp_limit_down = 0.01*resolution
                elif tech in ['IGCCAvgCF', 'Biopower']:
                    ramp_limit_up = 0.1*resolution
                    ramp_limit_down = 0.1*resolution
                elif tech in ['CCAvgCF','CTAvgCF']:
                    ramp_limit_up = min(0.6*resolution,1.0)
                    ramp_limit_down = min(0.6*resolution,1.0)
                else:
                    ramp_limit_up = 1.0
                    ramp_limit_down = 1.0
                
                # existing capacity
                if tech in generators.loc[bus].index:
                    p_nom = generators.at[bus, tech]
                else:
                    p_nom = 0.0

                if tech == 'LWR':
                    p_min_pu = 0.0
                    p_max_pu = 1.0  
                else:
                    p_max_pu = 1
                    p_min_pu = 0
                    
                extendable = tech in snakemake.config['extendable_techs']

                n.add(class_name="Generator",
                    name=f"{bus} {tech}",
                    bus=bus,
                    p_nom=p_nom,
                    p_nom_min=p_nom,
                    p_nom_extendable=extendable,
                    p_min_pu = p_min_pu,
                    p_max_pu = p_max_pu,
                    carrier=carrier,
                    capital_cost=item.capital_cost,
                    marginal_cost=item.marginal_cost,
                    lifetime=item.lifetime,
                    ramp_limit_down = ramp_limit_down,
                    ramp_limit_up = ramp_limit_up,
                    )
    return
              
              
def attach_storage(n, costs, generators):
    carriers = ["Batteries"]
    for bus in n.buses.index:
        for item in costs.itertuples():
            tech = item.techdetail
            carrier = item.technology_alias
            if carrier in carriers:
                # existing capacity
                if tech in generators.loc[bus].index:
                    p_nom = generators.at[bus, tech]
                else:
                    p_nom = 0.0
                    
                extendable = tech in snakemake.config['extendable_techs']

                n.add(class_name="StorageUnit",
                    name=f"{bus} {tech}",
                    bus=bus,
                    p_nom=p_nom,
                    p_nom_min=p_nom,
                    p_nom_extendable=extendable,
                    carrier=carrier,
                    capital_cost=item.capital_cost,
                    marginal_cost=item.marginal_cost,
                    lifetime=item.lifetime,
                    max_hours=float(tech.split(' ')[0].strip('Hr')),
                    cyclic_state_of_charge=False)
        
    return  


if __name__ == "__main__":
    
    n = pypsa.Network(snakemake.input.base_network)
    costs = load_costs()
    
    costs.to_csv("data/final_costs.csv")
    costs = costs.xs((slice(None), slice(None), model_years[0]))
    costs = costs.reset_index()
    costs.loc[costs['technology_alias']=='Solar', 'techdetail'] = 'Utility PV'
    
    generators = load_existing_generators()
    
    emissions = load_emissions()
    
    attach_load(n)
    add_carriers(n, 
                 costs=costs,
                 emissions=emissions)
    attach_renewables(n,
                      costs=costs,
                      generators=generators)
    attach_generators(n, 
                      costs=costs, 
                      generators=generators
                      )   
    attach_storage(n,
                   costs=costs,
                   generators=generators)
    
    # add co2 constraint
    n.add(class_name="GlobalConstraint",
          name="CO2 Limit",
          carrier_attribute='co2_emissions',
          sense="<=",
          investment_period=model_years[-1],
          constant=0) 
    
    n.export_to_netcdf(snakemake.output.elec_network)
    
    