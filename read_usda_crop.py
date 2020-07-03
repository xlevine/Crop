import sys, os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.figure_factory as ff
#import plotly.express as px
#from mpl_toolkits.basemap import Basemap,maskoceans

from os.path import expanduser
home = expanduser("~")

# data source: https://www.fsa.usda.gov/news-room/efoia/electronic-reading-room/frequently-requested-information/crop-acreage-data/index
# https://quickstats.nass.usda.gov/
# open and read xlsx file
def dic_usda(sheet_call):

    print(sheet_call)
    if sheet_call=='crop_acres_by_state':
        sheet_name = 'plant_and_fail'
    elif sheet_call=='failed_crop_acres_by_state':
        sheet_name = 'prevent'
    elif sheet_call=='total_acres_by_county':
        sheet_name = 'county_summary'
    elif sheet_call=='total_acres_by_crop':
        sheet_name = 'crop_summary'
    elif sheet_call=='crop_acres_by_county':
        sheet_name = 'county_data'
    elif sheet_call=='farms_by_county':
        sheet_name = 'all_farm_count'

    return sheet_name

def save_xlsx_usda_to_pickle(sheet_call):

    SHEET_NAME = dic_usda(sheet_call) 

    file_path = '/Users/xl332/Desktop/DataProjects/UsdaCrop/Data/'
    file_name = '2018_fsa_acres_120318.xlsx'
    file_df = pd.read_excel(file_path+file_name,sheet_name=SHEET_NAME)

    file_df.to_pickle(sheet_call + '.pkl')    

def read_var_usda(var_name, sheet_call):

    file_df = pd.read_pickle(sheet_call + '.pkl')
    header = list(file_df.columns.values)

    return file_df
    
def plot_var_usda(var_name, sheet_call):

    file_df = pd.read_pickle(sheet_call + '.pkl')
    fips = file_df['State County Code'].tolist()
    var = file_df[var_name].tolist()
        
    # remove rows with NaN in fips
    fips = fips[:-1]
    VAR = VAR[:-1]
    FIPS = [int(i) for i in fips]

    savefile_name = var_name + '_usa'
    fig = ff.create_choropleth(fips=FIPS, values=VAR)
    fig.layout.template = None
    fig.write_image(savefile_name+'.png')

def get_county_fips(state_code=0):

    file_df = pd.read_pickle('crop_acres_by_county' + '.pkl')

    if state_code==0:
        # 1. list all unique county codes
        fips = file_df['State County Code'].tolist()
    else:
        file_df_state = file_df.loc[file_df['State Code']==state_code]
        fips = file_df_state['State County Code'].tolist()

    fips_list = np.unique(fips)

    return fips_list

def read_crop_by_county(crop_name='All',Irrigation='All',Use='All'):

    # Output: Planted, Prevented, Failed Acres
    # Inputs: (1) All or specific crops
    #         (2) Irrigation practice: All, Yes, No
    #         (3)Intended use: All, or specific use

#    crop_list = ['WHEAT', 'OATS', 'CORN', 'SORGHUM']
    file_df = pd.read_pickle('crop_acres_by_county' + '.pkl')

    if crop_name!='All':
        crop_df = file_df.loc[file_df['Crop'] == crop_name]
    else:
        crop_df = file_df

    if Irrigation=='No':
        crop_df = crop_df.loc[file_df['Irrigation Practice'] == 'N']        
    elif Irrigation=='Yes':
        crop_df = crop_df.loc[file_df['Irrigation Practice'] == 'I']        

    if Use!='All':
        crop_df = crop_df.loc[file_df['Intended Use'] == Use]        

    # for every county, sum all contributions [other criteria can be added]
    crop_county_df = pd.DataFrame(columns = ['State','State County Code','Planted and Failed Acres'])

    FIPS = get_county_fips(0)

    # 2. iterates over each county code
#    index = 0
    for f in FIPS:
        crop_fips = crop_df.loc[crop_df['State County Code']==f]
        state_name_list = crop_fips['State'].to_list()
        if len(state_name_list)==0:
            state_name = ''
        else:
            state_name = state_name_list[0]
        crop_acre_all = np.sum(crop_fips['Planted and Failed Acres'])
        crop_county_n = pd.DataFrame([[state_name, int(f), crop_acre_all]],columns = ['State','State County Code','Planted and Failed Acres'])
        crop_county_df = pd.concat([crop_county_df,crop_county_n],ignore_index='True')
#        crop_county_df.loc[index] = [state_name, int(f), crop_acre_all]
#        index = index + 1

    return crop_county_df


def plot_crop_by_county(crop_name='All',Irrigation='All',Use='All'):

    crop_county_df = read_crop_by_county(crop_name,Irrigation,Use)
    fips = crop_county_df['State County Code'].tolist()
    acre = crop_county_df['Planted and Failed Acres'].tolist()

    savefile_name = crop_name + '_use_' + Use + '_irr_' + Irrigation + '_usa'
    fig = ff.create_choropleth(fips=fips, values=acre, scope=["usa"], title_text='Crop name')
    fig.layout.template = None
    fig.show()

#    fig.write_image(savefile_name+'.png')

#  goal = (1) load MODIS hdf data and plot it on continuous USA map (done with matplotlib).
#         (2) aggregate MODIS data at county-wide scale and plot it on county USA map (done with geopandas).
       
#  (a) Assign FIPS to each MODIS pixels: (i) read county geometry file: https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json , and (ii) check whether any pixel is located within a polygon (county) using this method: https://automating-gis-processes.github.io/CSC18/lessons/L4/point-in-polygon.html
#  (b) for every county, find corresponding pixels (and slices): (i) find most southward, northward, eastward and westward point of county polygon englobed within rectangle defined by extrema; (ii) find corresponding slice(s); (iii) within slice, locate pixels belonging to county polygon; create an array of slice & pixel indices for each county (allowing for fast reading and analysis).

#    granules for contiguous USA (500m SIN):
#    *h08v04*,*h09v04*,*h10v04*,*h11v04*,*h12v04*,*h13v04*
#    *h08v05*,*h09v05*,*h10v05*,*h11v05*,*h12v05*
#    *h08v06*,*h09v05*,*h10v05*

#    download above 14 granules for 1 year (e.g. 01/2018-12/2018) on dmr32 if size allows. Analysis may be done on PW9 if python modules exist. Start with yearly land use, and monthly/biweekly NDVI. Download more years depending on size requirements
