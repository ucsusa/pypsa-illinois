configfile: "config.yml"

from pathlib import Path
env_file = Path("./.env").resolve()
from dotenv import load_dotenv
load_dotenv(str(env_file))

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

rule build_topology:
    input: 
        supply_regions = "data/spatial_data/supply_regions.shp"
    output:
        buses = "data/buses.csv",   
        lines = "data/lines.csv"
    script: "scripts/build_topology.py"

# rule add_electricity:
#     input:
#         supply_regions = "data/spatial_data/supply_regions.shp",
#         load = "data/time_series/load.csv",
#         generators = "data/aggregated_generators.csv",
#         costs = "data/technology_costs.csv"
#     output:

