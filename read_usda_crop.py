import sys, os
from os import listdir
from os.path import isfile, join
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import plotly.figure_factory as ff

from os.path import expanduser
home = expanduser("~")

# data source: https://www.fsa.usda.gov/news-room/efoia/electronic-reading-room/frequently-requested-information/crop-acreage-data/index
# https://quickstats.nass.usda.gov/
data_file_path = '/Users/xl332/Desktop/DataProjects/UsdaCrop/Data/'
script_file_path = '/Users/xl332/Desktop/DataProjects/UsdaCrop/'

def find_filename(month,year):

    onlyfiles = [f for f in listdir(data_file_path) if isfile(join(data_file_path, f))]
    file_name_root = year + '_fsa_acres_' + month
    indices = [i for i, s in enumerate(onlyfiles) if file_name_root in s]
    file_name = onlyfiles[indices[0]]

    return file_name

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

sheet_call= 'crop_acres_by_county'
def save_xlsx_usda_to_pickle(year,month):

    SHEET_NAME = dic_usda(sheet_call) 
    file_name = find_filename(month,year)
    
    file_df = pd.read_excel(data_file_path+file_name,sheet_name=SHEET_NAME)

    file_df.to_pickle(sheet_call + '_' + year + month + '.pkl')    

def get_county_fips(year,month,state_code=0):

    file_df = pd.read_pickle(sheet_call + '_' + year + month + '.pkl')

    if state_code==0:
        # 1. list all unique county codes
        fips = file_df['State County Code'].tolist()
    else:
        file_df_state = file_df.loc[file_df['State Code']==state_code]
        fips = file_df_state['State County Code'].tolist()

    fips_list = np.unique(fips)

    return fips_list

def read_crop_by_county(year,month,crop_name='All',Irrigation='All',Use='All'):

    # Output: Planted, Prevented, Failed Acres
    # Inputs: (1) All or specific crops
    #         (2) Irrigation practice: All, Yes, No
    #         (3)Intended use: All, or specific use
    filename = sheet_call + '_' + year + month + '.pkl'
    onlyfiles = [f for f in listdir(script_file_path) if isfile(join(script_file_path, f))]
    if filename not in onlyfiles:
        save_xlsx_usda_to_pickle(year,month)

    file_df = pd.read_pickle(sheet_call + '_' + year + month + '.pkl')

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

    FIPS = get_county_fips(year,month,0)

    # 2. iterates over each county code
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

    # sum contributions from all counties
    total_planted = np.sum(crop_county_df['Planted and Failed Acres'])    
    crop_county_df['Planted and Failed Acres'] = crop_county_df['Planted and Failed Acres']/total_planted*100.0

    return crop_county_df

def plot_crop_by_county(year,month,crop_name='All',Irrigation='All',Use='All'):

    title_name = crop_name.lower()+ ' acreage in ' + month + '/' + year + ' (% of total acreage)'

    crop_county_df = read_crop_by_county(year,month,crop_name,Irrigation,Use)
    fips = crop_county_df['State County Code'].tolist()
    acre = crop_county_df['Planted and Failed Acres'].tolist()
    max_bar = round(np.amax(acre),1)
    nbar = 10
    binning_endpoints_list = list(np.arange(0,max_bar,max_bar/nbar))
    savefile_name = crop_name + '_use_' + Use + '_irr_' + Irrigation + '_usa'
    fig = ff.create_choropleth(fips=fips, values=acre, scope=["usa"], title_text=title_name, binning_endpoints = binning_endpoints_list)
    fig.layout.template = None
    fig.show()
