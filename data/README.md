# Data for PyPSA-Illinois

## Existing Power Plants

Data on existing power plants comes from EIA Form 823, accessed via the EIA APIv2 with the script `scripts/existing_power_plants.py`.


## Transmission Lines and Buses

The data for transmission lines and buses are available from the Energy Information Administration (EIA)
through the Homeland Infrastructure Foundation Level Database (HIFLD).

Unfortunately, due to the size of the file, the data are not included in this repository. Data about the
buses are imputed from the transmission line dataset. Users attempting to replicate this work should download
the data to the folder `pypsa-illinois/data/transmission-lines` and unzip. `pypsa-illinois` expects a shapefile dataset
called `pypsa-illinois/data/transmission_lines/Electric__Power_Transmission_Lines/Electric__Power_Transmission_Lines.shp`.

#### [Download the data here](https://atlas.eia.gov/datasets/geoplatform::transmission-lines/about)


## Electric Service Territories

The data for electric service territories is also availabe through EIA and HIFLD. Similar to the transmission line data,
there is no method for extracting the data programmatically. Users attempting to replicate this work should download
the data to the folder `pypsa-illinois/data/service-territories` and unzip. `pypsa-illinois` expects a shapefile dataset
called `pypsa-illinois/data/transmission_lines/Electric__Power_Transmission_Lines/Electric__Power_Transmission_Lines.shp`.

#### [Download the data here](https://atlas.eia.gov/datasets/geoplatform::electric-retail-service-territories-2/about)