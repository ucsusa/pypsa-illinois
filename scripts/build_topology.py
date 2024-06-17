# coding: utf-8

import networkx as nx
import pandas as pd
import geopandas as gpd
import numpy as np
from operator import attrgetter

def build_topology():
    regions = gpd.read_file(snakemake.input.supply_regions).set_index('name')
    centroids = regions\
                .to_crs(epsg=5070)\
                    .centroid\
                        .to_crs(epsg=4326)\
                            .to_frame()\
                                .rename(columns={0:'geometry'})
    
    line_config = snakemake.config['lines']
    v_nom = line_config['v_nom']

    buses = (regions.assign(v_nom=v_nom).drop('geometry', axis=1))

    # Lines from touching regions
    lines = pd.DataFrame([regions.index.values], columns=['bus0','bus1'])

    return buses, lines

if __name__ == "__main__":
    buses, lines = build_topology()

    buses.to_csv(snakemake.output.buses)
    lines.to_csv(snakemake.output.lines)
