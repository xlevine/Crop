This projects aims to map crops being grown in the USA based on (1) the USDA inventory of crop area, (2) a MODIS vegetation map. 

# Mapping Land Surface using MODIS datasets  

I.  First, we map MODIS vegetation for any given geographical unit (County or State). To do so, you need to follow the following recipe:

1. download read_modis.py and read_hdf_file.py in this code respository to your local working directory

2. download MODIS tiles for a dataset of interest (see instructions at the end of this tutorial). 

3. fix all hardcoded paths to fit your local directory tree structure. 

4. open a Python3 prompt

5. in the prompt: "import read_modis" (this script load when all required packages are installed: matplotlib, numpy, pandas, shapely, json, cartopy, pyhdf.SD)

6. in the prompt: "from read_modis import map_landuse_2D"

7. in the prompt: "map_landuse_2D(year,location)" . Both inputs are strings (e.g. year='2017'); the location string either refers to the FIPS code of the county (e.g '06037' for Los Angeles county) or the State ID (e.g. 'CA' for California).  

FIPS code may be found here: https://www.nrcs.usda.gov/wps/portal/nrcs/detail/national/home/?cid=nrcs143_013697

After completing step (7), you will see a PNG file appear in your local directory, named after the fips code or State ID (e.g. CA.png or 06037.png). 

By default, the code generate plots for the yearly land use ('LC_Type_1') based on the Annual International Geosphere-Biosphere Programme (IGBP) classification. This classification bins land use into 17 classes of 'natural' and 'anthropogenic' landscapes.


<img src="https://github.com/xlevine/Crop/blob/master/plots/CA.png" width="300"><img src="https://github.com/xlevine/Crop/blob/master/plots/IL.png" width="300"><img src="https://github.com/xlevine/Crop/blob/master/plots/CT.png" width="300">

<img src="https://github.com/xlevine/Crop/blob/master/plots/06037.png" width="300"><img src="https://github.com/xlevine/Crop/blob/master/plots/17031.png" width="300"><img src="https://github.com/xlevine/Crop/blob/master/plots/09009.png" width="300">

To download MODIS data: 

1. go to https://lpdaac.usgs.gov/product_search/

2. search for "MODIS"

3. select dataset: e.g. "MCD12Q1 v006" for Land Cover Type Yearly, "MCD12Q2 v006" for Land Cover Dynamics Yearly; "MCD15A2H v006" for Leaf Area Index/FPAR 8-Day.

4. determine variable of interest (e.g. 'LC_Type_1' for land use when using the MCD12Q1 dataset). 

5. Follow the download link to the Nasa Earthdata Search. Select all revelevant files by year (e.g. 2017) and tiles. MODIS data is outputed in 3 grid format: Climate Modelling (CMG), Lambert Azimuthal Equal-Area Tile (Polar), and Sinusoidal Tile (SIN) grids. 

For more information about the different grids:  https://modis-land.gsfc.nasa.gov/MODLAND_grid.html 

This code can handle data in SIN grids at all available resolution (100m, 250m, 500m, 1000m), which is the most common format for high-resolution datasets. This code can handle any number of temporal resolutions (yearly, monthly, 14 days, 8 days, 4 days, daily).

For data in SIN grid, select tiles of interest (e.g. h9-h12, v4-v6) and download them using the wget file provided by the Earthdata download tool.

Note: you will need to update the variable name (in read_modis.py), and the dataset name and variable list if not already predefined (in dic_modis.py). To produce a map of your chosen variable, simply input the time frame ('YYYY' for a yearly dataset, 'YYYYMM' for a monthly-mean one, 'YYYYMMDD' for a submonthly one) and the location ID (county / state) in map_landuse_2D(). Et voila!

# Mapping Crop acreage using USDA datasets  

https://www.fsa.usda.gov/news-room/efoia/electronic-reading-room/frequently-requested-information/crop-acreage-data/index 

read_usda_crop.py

