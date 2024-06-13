configfile: "config.yml"

from pathlib import Path
env_file = Path("./.env").resolve()

rule retrieve_supply_regions:
    output: "data/spatial_data/supply_regions.shp"
    script: "scripts/retrieve_supply_regions.py"