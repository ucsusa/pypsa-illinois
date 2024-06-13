import geopandas as gpd


if __name__ == "__main__":

    state = snakemake.config['state_abbr']

    url = (f"https://services1.arcgis.com/Hp6G80Pky0om7QvQ/arcgis/"
           f"rest/services/Retail_Service_Territories/FeatureServer/"
           f"0/query?where=STATE%20%3D%20'{state}'&outFields=STATE,"
           f"CNTRL_AREA,PLAN_AREA,HOLDING_CO&outSR=4326&f=json")

    gdf = gpd.read_file(url)
    gdf = gdf.dissolve('CNTRL_AREA').reset_index(drop=False)
    
    if state == 'IL':
        gdf = gdf.loc[[1,2]].reset_index(drop=True)
        
        gdf = gpd.GeoDataFrame(geometry=[gdf.loc[0].geometry.difference(gdf.loc[1].geometry), 
                                         gdf.loc[1].geometry]).set_crs(epsg=4326)
        gdf['name'] = ['MISO-Z4','ComEd']
        gdf = gdf.set_index('name')
        
    gdf.to_file("data/spatial_data/supply_regions.shp")
        