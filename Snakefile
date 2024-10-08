configfile: "config.yml"

from pathlib import Path
env_file = Path("./.env").resolve()
from dotenv import load_dotenv
load_dotenv(str(env_file))

if config['myopic']:
    results_folder = f"{config['scenario']}_v{config['version']}_myopic"
else:
    results_folder = f"{config['scenario']}_v{config['version']}"

rule targets:
    input:
        f"results/{results_folder}/figures/illinois_dispatch.png",
        f"results/{results_folder}/networks/illinois_solved.nc",
        f"dag.png"
    # input:
    #     expand("results/")

rule retrieve_supply_regions:
    output: 
        supply_regions = "data/spatial_data/supply_regions.shp"
    script: "scripts/retrieve_supply_regions.py"

rule retrieve_costs:
    output: 
        costs = "data/technology_costs.csv",
        heatrates = "data/heatrates.csv"
    script: "scripts/retrieve_costs.py"

rule retrieve_fuel_costs:
    input:
        heatrates = "data/heatrates.csv"
    output: 
        fuel_costs = "data/thermal_fuel_costs.csv",
        fuel_cost_timeseries = "data/thermal_fuel_cost_ts.csv"
    script: "scripts/retrieve_fuel_prices.py"

rule retrieve_load:
    output:
        load = "data/time_series/load.csv"
    script: "scripts/retrieve_load.py"

rule retrieve_existing_generators:
    output: 
        generators = "data/existing_generators.csv",
        aggregated = "data/aggregated_generators.csv",
        build_year = "data/build_year_generators.csv"
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
        config_file = 'config.yml',
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
        fuel_cost_timeseries = "data/thermal_fuel_cost_ts.csv",
        wind_profile = "data/time_series/wind.csv",
        solar_profile = "data/time_series/solar.csv",
        base_network = "data/networks/base_network.nc",
        emissions = "data/technology_emissions.csv",
        build_years = "data/build_year_generators.csv"
    output:
        elec_network = "data/networks/electricity_network.nc",
        final_costs = "data/final_costs.csv"
    script: "scripts/add_electricity.py"

rule solve_network:
    input:
        elec_network = "data/networks/electricity_network.nc"
    output:
        solved_network = f"results/{results_folder}/networks/illinois_solved.nc"
    script: "scripts/solve_network.py"

rule plot_results:
    input:
        solved_network = f"results/{results_folder}/networks/illinois_solved.nc"
    output: 
        dispatch_figure = f"results/{results_folder}/figures/illinois_dispatch.png",
        capacity_figure = f"results/{results_folder}/figures/illinois_capacity.png",
        emissions_figure = f"results/{results_folder}/figures/illinois_emissions.png",
        active_units_figure = f"results/{results_folder}/figures/illinois_active_units.png",
        monthly_generation_figure = f"results/{results_folder}/figures/illinois_monthly_generation.png"
    script: "scripts/plot_results.py"

rule build_dag:
    input: "Snakefile"
    output:
        "dag.png"
    shell:
        "snakemake --dag | dot -Tpng > {output}"