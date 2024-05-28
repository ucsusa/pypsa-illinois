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


## Solar and Wind Resource Data 

Solar and wind resource data are extracted from the National Solar Radiation Database (NSRDB) and Wind Toolkit (WTK), respectively.

The data are extracted using functions in the file `functions/nrel_data_api.py`. These functions are deployed in the notebook `notebooks/nrel_data_download.ipynb`.

### [NSRDB](https://nsrdb.nrel.gov/about/what-is-the-nsrdb)

> The NSRDB (Sengupta et al., 2018) is a high temporal and spatial resolution dataset consisting of the three most widely used measurements of solar radiation—global horizontal, direct normal, and diffuse horizontal irradiance—as well as other meteorological data. The earlier versions of the NSRDB were modeled using cloud and weather information primarily collected at airports. An adequate number of locations and temporal and spatial scales were used to accurately represent regional solar radiation climates. More details on the NSRDB version history can be found here. The current NSRDB is modeled using the NREL’s Physical Solar Model (PSM) with inputs from multi-channel measurements obtained from the Geostationary Operational Environmental Satellite (GOES) of the National Oceanic and Atmospheric Administration (NOAA), the Interactive Multisensor Snow and Ice Mapping System (IMS) of the National Ice Center (NIC), and the Moderate Resolution Imaging Spectroradiometer (MODIS) and Modern Era Retrospective analysis for Research and Applications, version 2 (MERRA-2), of the National Aeronautics and Space Administration (NASA). The PSM is a two-step physical modeling process, in which cloud and aerosol properties are derived, collected and resampled in the initial step and then fed as inputs into a radiative transfer model, the Fast All-sky Radiation Model for Solar applications (FARMS) (Xie et al., 2016), in the subsequent step. Using the FARMS with Narrowband Irradiances for Tilted surfaces (FARMS-NIT) (Xie and Sengupta, 2018; Xie et al., 2018; Xie et al., 2019), the NSRDB can also provide users with spectral-on-demand irradiances based on their selection of time, location, and PV orientation.


### Wind Toolkit

> The Wind Integration National Dataset (WIND) Toolkit is an update and expansion of the Eastern Integration Data Set and Western Wind Integration Data Set. It supports the next generation of wind integration studies. The WIND Toolkit includes meteorological conditions and turbine power for more than 126,000 sites in the continental United States for the years 2007—2013. Read more at http://www.nrel.gov/grid/wind-toolkit.html