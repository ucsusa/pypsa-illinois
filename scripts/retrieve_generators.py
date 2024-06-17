import numpy as np
import pandas as pd
import requests
import io
import os
import sys
import json
import geopandas as gpd
from gridstatus import EIA
from dotenv import load_dotenv
from pathlib import Path


curr_dir_os = Path(os.path.dirname(os.path.abspath(__file__)))

load_dotenv(curr_dir_os/"../.env")


eia = EIA()

BASE_URL = "https://api.eia.gov/v2/"

dataset = "electricity/operating-generator-capacity"
route = dataset + "/data"


DATASET_CONFIG = {
    "electricity/operating-generator-capacity": {
        "index": [
            "period",
            "balancing_authority_code",
            "stateid",
            "plantid",
            # "balancing_authority_code",
        ],
        # "handler": _handle_rto_interchange,
    },
}

start_str = "2024-01"
end_str = "2024-04"
frequency= 'monthly'
facets = {
        #   'balancing_authority_code':['SWPP'],
          'stateid':['IL'],
          'status':['OP','SB','OA']
          }

params = {
            "start": start_str,
            "end": end_str,
            "frequency": frequency,
            "facets": facets,
            "offset": 0,
            "length": 5000,
            "data":["county", 
                    "longitude",
                    "latitude",
                    "nameplate-capacity-mw",
                    "operating-year-month",
                    "planned-retirement-year-month"],
            # pagination breaks if not sorted because
            # api doesn't return in stable order across requests
            "sort": [
                {"column": col, "direction": "asc"}
                for col in DATASET_CONFIG[dataset]["index"]
            ],
        }
headers = {
            "X-Api-Key": eia.api_key,
            "X-Params": json.dumps(params),
        }


if __name__ == "__main__":
    
    r = requests.get(BASE_URL+route, headers=headers)
    
    if r.status_code == 200:
        print("Successfully downloaded data from EIA")
    else:
        raise requests.HTTPError()
    
    response = r.json()['response']
    
    print(f"Total records: {response['total']}")
    
    data = response['data']
    
    df = pd.DataFrame(data)
    df.drop(columns=['stateName','sector','sectorName','entityid',
                     'balancing-authority-name','statusDescription', 
                     'entityName','nameplate-capacity-mw-units',
                     'energy_source_code','generatorid','status',
                     'unit'
                     ],
            inplace=True)
    
    locations = gpd.points_from_xy(x=df['longitude'], y=df['latitude'], crs='EPSG:4326')
    
    gdf = gpd.GeoDataFrame(df, geometry=locations)
    
    gdf['nameplate-capacity-mw'] = gdf['nameplate-capacity-mw'].astype('float')
    
    gdf['balancing_authority_code'] = gdf['balancing_authority_code'].replace({"MISO":"MISO-Z4",'PJM':"ComEd"})
    
    gdf.to_csv(snakemake.output.generators)
    
    scale = snakemake.config['geo_res']
    idx_opts = {"rto":"balancing_authority_code",
                "county":"county"}
    
    gen_agg = gdf.pivot_table(index=idx_opts[scale],
                                columns='technology',
                                values='nameplate-capacity-mw',
                                aggfunc='sum')
    
    gen_agg.to_csv(snakemake.output.aggregated)