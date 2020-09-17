# Crop 

This projects aims to map crops being grown in the USA based on (1) the USDA inventory of crop area, (2) a MODIS vegetation map. 

To plot land cover type in a county, follow the following recipe:

open a Python3 prompt

import read_modis_to_county (this script load when all required packages are installed)

from read_modis_to_county import plot_data_county_2D

plot_data_county_2D('LC_Type1',year,fips_county)

where year and fips_county are strings.
