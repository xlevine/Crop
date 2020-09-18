# Mapping Land Surface using MODIS datasets  

This projects aims to map crops being grown in the USA based on (1) the USDA inventory of crop area, (2) a MODIS vegetation map. 

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

5. update read_modis.py and read_hdf_file.py to update the variable name, the dataset name, and the spatial resolution if necessary. To plot this new quantity, simply input the time frame ('YYYY' for a yearly dataset, 'YYYYMM' for a monthly-mean one, 'YYYYMMDD' for a submonthly one) and the location ID (county / state) in map_landuse_2D(). Et voila!
