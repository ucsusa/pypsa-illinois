import networkx as nx
import pandas as pd
import numpy as np
import scipy as sp
import geopandas as gpd
from tqdm import tqdm
import pypsa

model_year = snakemake.config['model_year']

def annuity(r, n):
    return r / (1-1/(1+r)**n)


def create_snapshots():
    snapshots = pd.date_range(start=f"{model_year}-01-01", 
                            end=f"{model_year+1}-01-01", 
                            freq='1h', 
                            inclusive='left')
    
    return snapshots


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
    
    load = load.loc[str(snakemake.config['load_year'])]
    
    snapshots = create_snapshots()
    
    load.set_index(snapshots, inplace=True)
    
    load = load.resample(f"{snakemake.config['time_res']}h").sum()
    
    print('Adding loads to model')
    for bus in tqdm(n.buses.index):
        n.add(class_name="Load",
              name=bus,
              bus=bus,
              p_set=load[bus]
              )
    
    return


def attach_generators(n, costs, generators):
    costs = costs.xs((slice(None), slice(None), model_year))

    available_carriers = costs.index.get_level_values('technology_alias').unique().to_list()
    
    snapshots = create_snapshots()
    wind_profile = pd.read_csv(snakemake.input.wind_profile, parse_dates=True, index_col=0)
    wind_profile.set_index(snapshots, inplace=True)
    wind_profile = wind_profile.resample(f"{snakemake.config['time_res']}h").mean()
    
    solar_profile = pd.read_csv(snakemake.input.solar_profile, parse_dates=True, index_col=0)
    solar_profile.set_index(snapshots, inplace=True)
    solar_profile = solar_profile.resample(f"{snakemake.config['time_res']}h").mean()
    
    # add carriers
    for carrier in available_carriers:
        n.add(class_name="Carrier",
              name=carrier, 
              color=snakemake.config['carrier_colors'][carrier])
        
    # flatten index
    costs = costs.reset_index()
    costs.loc[costs['technology_alias']=='Solar', 'techdetail'] = 'Utility PV'
    for bus in n.buses.index:
        for item in costs.itertuples():
            tech = item.techdetail
            carrier = item.technology_alias
            capital_cost = item.capital_cost
            marginal_cost = item.marginal_cost
            lifetime = item.lifetime
            
            if carrier == 'Batteries':
                class_name = "StorageUnit"
                max_hours = float(tech.split(' ')[0].strip('Hr'))
                cyclic_state_of_charge=True
            else: 
                class_name = "Generator"
                max_hours = 0.0
                cyclic_state_of_charge=False
            
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
            else:
                p_max_pu = 1
                
            extendable = tech in snakemake.config['extendable_techs']

            n.add(class_name=class_name,
                  name=f"{bus} {tech}",
                  bus=bus,
                  p_nom=p_nom,
                  p_nom_min=p_nom,
                  p_max_pu=p_max_pu,
                  p_nom_extendable=extendable,
                  carrier=carrier,
                  capital_cost=capital_cost,
                  marginal_cost=marginal_cost,
                  lifetime=lifetime,
                  max_hours=max_hours,
                  cyclic_state_of_charge=cyclic_state_of_charge)
    return
                


if __name__ == "__main__":
    
    n = pypsa.Network(snakemake.input.base_network)
    costs = load_costs()
    
    costs.to_csv("data/final_costs.csv")
    
    generators = load_existing_generators()
    
    attach_load(n)
    attach_generators(n, costs=costs, generators=generators)    
    
    n.export_to_netcdf(snakemake.output.elec_network)
    
    