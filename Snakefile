configfile: "config.yml"

from pathlib import Path
env_file = Path("./.env").resolve()
from dotenv import load_dotenv
load_dotenv(str(env_file))

rule targets:
    input:
        "results/figures/illinois_dispatch.png",
        "results/networks/illinois_solved.nc"

rule retrieve_supply_regions:
    output: 
        supply_regions = "data/spatial_data/supply_regions.shp"
    script: "scripts/retrieve_supply_regions.py"

rule retrieve_costs:
    output: 
        costs = "data/technology_costs.csv"
    script: "scripts/retrieve_costs.py"

rule retrieve_load:
    output:
        load = "data/time_series/load.csv"
    script: "scripts/retrieve_load.py"

rule retrieve_existing_generators:
    output: 
        generators = "data/existing_generators.csv",
        aggregated = "data/aggregated_generators.csv"
    script: "scripts/retrieve_generators.py"

rule retrieve_renewable_profiles:
    input:
        supply_regions = "data/spatial_data/supply_regions.shp"
    output:
        wind = "data/time_series/wind.csv",
        solar = "data/time_series/solar.csv"
    script: "scripts/retrieve_renewables.py"

rule retrieve_co2_emissions:
    output: 
        emissions = "data/technology_emissions.csv"
    script: "scripts/retrieve_emissions.py"

rule build_topology:
    input: 
        supply_regions = "data/spatial_data/supply_regions.shp"
    output:
        buses = "data/buses.csv",   
        lines = "data/lines.csv"
    script: "scripts/build_topology.py"

rule build_base_network:
    input: 
        buses='data/buses.csv',
        lines='data/lines.csv'
    output: 
        network="data/networks/base_network.nc"
    script: "scripts/build_base_network.py"

rule add_electricity:
    input:
        supply_regions = "data/spatial_data/supply_regions.shp",
        load = "data/time_series/load.csv",
        generators = "data/aggregated_generators.csv",
        costs = "data/technology_costs.csv",
        wind_profile = "data/time_series/wind.csv",
        solar_profile = "data/time_series/solar.csv",
        base_network = "data/networks/base_network.nc",
        emissions = "data/technology_emissions.csv"
    output:
        elec_network = "data/networks/electricity_network.nc",
        final_costs = "data/final_costs.csv"
    script: "scripts/add_electricity.py"

rule solve_network:
    input:
        elec_network = "data/networks/electricity_network.nc"
    output:
        solved_network = "results/networks/illinois_solved.nc"
    script: "scripts/solve_network.py"

rule plot_results:
    input:
        solved_network = "results/networks/illinois_solved.nc"
    output: 
        dispatch_figure = "results/figures/illinois_dispatch.png",
        capacity_figure = "results/figures/illinois_capacity.png",
        emissions_figure = "results/figures/illinois_emissions.png"
    script: "scripts/plot_results.py"