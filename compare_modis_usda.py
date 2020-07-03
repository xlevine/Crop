import sys, os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy import stats
import plotly.figure_factory as ff
import json
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
from read_modis_to_county import read_modis_landuse
from read_usda_crop import read_crop_by_county, save_xlsx_usda_to_pickle
from os.path import expanduser
home = expanduser("~")

# you have spent enough time on this regression problem. I want to move on and start producing: 
# (a) main crops distribution by county
# (b) land use by pixel, county, state
# (c) decent match in cropland area between MODIS and USDA

# ("crop regime", to group counties by production regime): this is doable.
# (1) list main crops in each county from USDA, and compare amount explained when comparing with MODIS. [eqv. "regression"]
# (2) cluster cropland pixels according to their NDVI timeseries (one or more years); compute size of clusters and tie them to USDA crop type.

#np.set_printoptions(threshold=sys.maxsize)

state_list = ["Alabama", "Arkansas", "Arizona", "California", "Colorado", "Connecticut","Delaware", "Florida", "Georgia", "Iowa", "Idaho", "Illinois", "Indiana", "Kansas", "Kentucky", "Louisiana", "Massachusetts", "Maryland", "Maine", "Michigan", "Minnesota", "Missouri", "Mississippi", "Montana", "North Carolina", "North Dakota", "Nebraska", "New Hampshire", "New Jersey", "New Mexico", "Nevada", "New York", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Virginia", "Vermont", "Washington", "Wisconsin", "West Virginia", "Wyoming"]
tickslabel = ['Evgrn Ndlf','Evgrn Brdlf','Decdu Ndlf','Decdu Brdlf','Mixd Forest','Closd Shrub','Opn Shrub','Woody Svnas','Svnas','Grasslnds','Wetlnds','Crop','Urban','Crop/Natrl','Snow/Ice','Barren','Water']
#crop_list = ['BARLEY','BEANS','CANOLA','CHICKPEAS','CORN','COTTON','FLAXSEED','HAY','HOPS','LENTILS','MILLET','MINT','MUSTARD','OATS','PEANUTS','PEAS','RAPESEED','RICE','RYE','SAFFLOWER','SORGHUM','SOYBEANS','SUGARBEETS','SUGARCANE','SUNFLOWER','TARO','TOBACCO','WHEAT']
crop_list = ['ALFALFA','BARLEY','BEANS','CANOLA','CORN','COTTON UPLAND','COVER CROP','LENTILS','MILLET','OATS','PEANUTS','POTATOES','RICE','RYE','SORGHUM','SORGHUM FORAGE','SORGHUM DUAL PURPOSE','SOYBEANS','SUGARBEETS','SUGARCANE','SUNFLOWER','TRITICALE','WHEAT']

Vegetable_list = ['ARTICHOKES','ASPARAGUS','BEANS','BEETS','BROCCOLI','BRUSSELS SPROUTS','CABBAGE','CARROTS','CAULIFLOWER','CELERY','CUCUMBERS','EGGPLANT', 'GARLIC','GINGER ROOT','GREENS','LETTUCE','MELONS','OKRA','ONIONS','PEAS','PEPPERS','PICKLES','POTATOES','PUMPKINS','RADISHES','SPINACH','SQUASH','SWEET CORN','SWEET POTATOES','TOMATOES']

# determine which crops cover may be consistent with MODIS estimated crop area
# 1. regress total with each predictor, compute R-2
# 2. rank predictor by their R-2
# 3. construct model by adding predictor, starting by that with the greatest R-2, and discarding those that do not improve the adjusted R-2 of the multilinear regression (treating every crop as independent variables)
# 4. Potential Issue: when dealing with crops only important to certain areas, I risk discarding them because they do not improve the adj R-2 (?); this could be remediated by performing regressions over important agricultural areas (as determined by the USDA database): perform steps above; all variables deemed significant for each region will be kept for all.   
# Note on #4: when determining important agricultural areas, we could either use all counties in States or we could aggregate contiguous counties sharing similar productive profiles (K-mean clusters?).


def rank_crop(state_name):

    # read full list of crop (excel table)
    sheet_call = 'total_acres_by_crop'
#    save_xlsx_usda_to_pickle(sheet_call)
    file_df = pd.read_pickle(sheet_call + '.pkl')
    header = list(file_df.columns.values)
#    print(header)
    var = file_df['State'].tolist()
    crop_list = var[3:-2]
#    acre = file_df['(All)'].tolist()
#    acre_list = acre[3:-2]
#    k_list = np.argsort(-np.asarray(acre_list))
##    print(np.array(crop_list)[k_list])
#    print(np.array(acre_list)[k_list])
#    crop_array = np.array(crop_list)[k_list]
#    acre_array = np.array(acre_list)[k_list]
#    crop_array = crop_array[np.where(acre_array>=1.0e5)]
#    acre_array = acre_array[np.where(acre_array>=1.0e5)]

    crop_list = ['SOYBEANS', 'CORN', 'WHEAT', 'COTTON  UPLAND', 'ALFALFA', 'RICE', 'BARLEY', 'BEANS', 'OATS', 'CANOLA', 'PEANUTS', 'SUNFLOWERS', 'SUGAR BEETS', 'PEAS', 'MILLET', 'TRITICALE', 'SUGARCANE', 'POTATOES', 'RYE', 'LENTILS', 'ALMONDS', 'PECANS', 'COTTON  ELS', 'TOMATOES', 'FLAX', 'TOBACCO FLUE CURED', 'CLOVER', 'ORANGES', 'SAFFLOWER', 'PISTACHIOS', 'POTATOES SWEET', 'GRAPES', 'SESAME', 'ONIONS', 'WALNUTS']

#    print(k_list)
#    crop_list = ['ALFALFA','BARLEY','BEANS'] 

# for each item in the list, regress with total MODIS; store r-2
#    compute_usda = 'False'

# diagnose crop area for each crop in each county from USDA
#    compute_usda = 'True'
    compute_usda = 'False'
    if compute_usda=='True':
        acreage = []
        for crop_name in crop_list:
            print(crop_name)
            crop_county_df = read_crop_by_county(crop_name)
            if state_name!='USA':
                crop_state_df = crop_county_df.loc[crop_county_df['State']==state_name]
                acre = crop_state_df['Planted and Failed Acres'].tolist()
            else:
                acre = crop_county_df['Planted and Failed Acres'].tolist()
            acreage.append(acre)
        np.save('acreage_' + state_name + '.npy',[crop_list,acreage]) 
    else:
        [crop_list,acreage] = np.load('acreage_' + state_name + '.npy',allow_pickle=True)  


    # diagnose crop area in each county from MODIS
    landuse_county_df = load_modis_landuse(state_name)
    acre_luse_ = landuse_county_df['Crop'] + landuse_county_df['Crop/Natrl']
    fips_luse = landuse_county_df['fips'].tolist()

    # get fips in USDA dataset (could be merge with above loop)
    crop_county_df = read_crop_by_county()
    if state_name!='USA':
        crop_state_df = crop_county_df.loc[crop_county_df['State']==state_name]
        fips_usda = crop_state_df['State County Code'].tolist()
    else:
        fips_usda = crop_county_df['State County Code'].tolist()

    acreage_modis = []
    for fips in fips_usda:
        k_fips = np.argwhere(np.asarray(fips_luse)==fips)
#        acre_usda = acreage[n] 
        if len(k_fips)!=0:
            acre_luse = acre_luse_[k_fips[0][0]]
        else:
            acre_luse = np.nan
        acreage_modis.append(acre_luse)

    acreage_modis = np.asarray(acreage_modis)
    acreage_modis[np.where(acreage_modis==0)] = np.nan
    k_modis = np.squeeze(np.argwhere(np.isnan(acreage_modis)==False))

    acreage = np.asarray(acreage)
    r_squared = []; p_values = []
    for N in np.arange(0,len(crop_list),1):
        acreage_usda = np.squeeze(acreage[N])
        acreage_usda[np.where(acreage_usda==0)] = np.nan
#        print(acreage_usda)
        k_usda = np.squeeze(np.argwhere(np.isnan(acreage_usda)==False))
        k_crop = np.intersect1d(k_usda,k_modis)
        acre_usda = acreage_usda[k_crop]
        acre_modis = acreage_modis[k_crop]
        if len(acre_modis)==0:
            slope=np.nan; intercept=np.nan; r_value=np.nan; p_value=np.nan; std_err=np.nan
        else:
            slope, intercept, r_value, p_value, std_err = stats.linregress(acre_usda,acre_modis)
#        print(crop_list[N],'R-2='+str(round(r_value**2,2)), ' with sample size ', len(k_crop))
        if len(k_crop)<=50:
            slope=np.nan; intercept=np.nan; r_value=np.nan; p_value=np.nan; std_err=np.nan            
        r_squared.append(r_value**2)
        p_values.append(p_value)
        
    p_values = np.asarray(p_values)
    r_squared = np.asarray(r_squared)

#    kr = np.argsort(-r_squared)
    kr = np.argsort(p_values)

    r_squared_sorted = r_squared[kr]
    crop_list_sorted = crop_list[kr]
    p_values_sorted = p_values[kr]

    r_squared_sorted[np.where(p_values_sorted>=0.05)] = np.nan
    p_values_sorted[np.where(p_values_sorted>=0.05)] = np.nan
    crop_list_sorted[np.where(p_values_sorted>=0.05)] = np.nan

    r_squared_sorted = r_squared_sorted[np.where(np.isnan(r_squared_sorted)==False)]
    crop_list_sorted = crop_list_sorted[np.where(np.isnan(r_squared_sorted)==False)]
    p_values_sorted = p_values_sorted[np.where(np.isnan(r_squared_sorted)==False)]

    # here we build a model as follows: start with first crop, perform regression, compute R-2
    acreage_usda = np.nan * np.zeros(np.shape(acreage[0]))
    r2_adj = 0
    crop_list_ = []
    for N in np.arange(0,len(crop_list_sorted),1):
        acreage_usda[np.isnan(acreage_usda)==True] = 0
        acreage_usda_0 = acreage_usda
        r2_adj_0 = r2_adj 
        acreage_usda_new = np.squeeze(acreage[kr[N]])
        acreage_usda = np.add(acreage_usda,acreage_usda_new)
        acreage_usda[np.where(acreage_usda==0)] = np.nan
        k_usda = np.squeeze(np.argwhere(np.isnan(acreage_usda)==False))
        k_crop = np.intersect1d(k_usda,k_modis)
        acre_usda = acreage_usda[k_crop]
        acre_modis = acreage_modis[k_crop]
        slope, intercept, r_value, p_value, std_err = stats.linregress(acre_usda,acre_modis)
        r2_value = np.power(r_value,2)
        r2_adj = 1 - ((1-r2_value)*(len(k_crop)-1)/(len(k_crop)-len(crop_list_)-1))
        print(crop_list_sorted[N],'R-2='+str(r2_adj), ' with sample size ', len(k_crop))
        if r2_adj < r2_adj_0:
            acreage_usda = acreage_usda_0
            r2_adj = r2_adj_0
            print(crop_list_sorted[N] + ' is discarded')
        else:
            crop_list_.append(crop_list_sorted[N])
            print(crop_list_sorted[N] + ' is retained')

    print(crop_list_sorted)
    print(crop_list_)
            
    return crop_list_sorted

    # rank variable from their R-2 (discard all insignificant R-2)

# maybe: for each county list crops amounting for at least 10% of the total amount of cropland (excluding grass, mixed forage, cover crop, garden home, sorghum forage, sorghum, skip rows, rye, tree timber, etc)

def load_usda_crop(state_name,crop_name):

    if crop_name=='All':
        crop_county_df = read_crop_by_county(crop_name='All',Irrigation='All',Use='All')    
    else:
        crop_county_df = read_crop_by_county(crop_name,Irrigation='All',Use='All')

    if state_name!= 'USA':            
        crop_county_df = crop_county_df.loc[crop_county_df['State']==state_name]

    crop_acreage_df = pd.DataFrame(columns = ['fips', 'acres'])
    acre_crop_ = crop_county_df['Planted and Failed Acres'] #.tolist()
    fips_crop_ = crop_county_df['State County Code'] #.tolist()        
    fips_list = []
    for fips in fips_crop_:
        acreage = np.float(acre_crop_.loc[fips_crop_==fips].to_string(index=False))
        acreage = acreage * acre_to_ha
        crop_acreage = pd.DataFrame([[fips,acreage]],columns = ['fips', 'acres'])
        if fips not in fips_list:
            crop_acreage_df = pd.concat([crop_acreage_df,crop_acreage],ignore_index='True')
            fips_list.append(fips) 
        else:
            acreage_new = crop_acreage_df['acres'].loc[crop_acreage_df['fips']==fips] + acreage
            crop_acreage_df['acres'].loc[crop_acreage_df['fips']==fips] = acreage_new

    return crop_acreage_df

acre_to_ha  = 0.40468564224

def load_modis_landuse(state_name):
    
    if state_name=='USA':
        landuse_county_df = pd.DataFrame(columns = np.concatenate((['fips'],tickslabel),axis=0))
        for state_name in state_list:
            landuse_county_df_state = read_modis_landuse(state_name,'LC_Type1',timestamp='2017')
            landuse_county_df = pd.concat([landuse_county_df,landuse_county_df_state],ignore_index='True')
    else:
        landuse_county_df = read_modis_landuse(state_name,'LC_Type1',timestamp='2017')

    return landuse_county_df

def compare_acres_modis_usda_by_county(state_name='USA'):
    
    # time stamp info must be added in USDA file
    crop_acreage_df = load_usda_crop(state_name)
    landuse_county_df = load_modis_landuse(state_name)
    acre_crop_ = crop_acreage_df['acres'].tolist()
    fips_crop_ = crop_acreage_df['fips'].tolist()

    acre_luse = landuse_county_df['Crop'] + landuse_county_df['Crop/Natrl']
#    acre_luse = landuse_county_df['Crop']# + landuse_county_df['Crop/Natrl']
    fips_luse = landuse_county_df['fips']

    fips_crop_ = [int(i) for i in fips_crop_] 
    farmed_acreage = pd.DataFrame(columns = ['fips', 'USDA', 'MODIS'])

    n = 0
    for fips in fips_luse:
        k_fips = np.argwhere(np.asarray(fips_crop_)==fips)
        acre_modis = acre_luse[n] 
        if len(k_fips)!=0:
            acre_usda = acre_crop_[k_fips[0][0]]
        else:
            acre_usda = np.nan
        farmed_acreage.loc[n] = [int(fips),acre_usda,acre_modis]
        n = n + 1

    return farmed_acreage

def regress_modis_usda(state_name):

    farmed_acreage = compare_acres_modis_usda_by_county(state_name)

    x = farmed_acreage['USDA']
    y = farmed_acreage['MODIS']
    x = np.asarray(x)
    y = np.asarray(y)
    kx_p = np.squeeze(np.argwhere(np.isnan(x)==False))
    ky_p = np.squeeze(np.argwhere(np.isnan(y)==False))
    k_p = np.intersect1d(kx_p,ky_p)
    x_p = x[k_p]
    y_p = y[k_p]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_p,y_p)
    print(r_value**2)

    fig = plt.figure()
    plt.plot(x,y,'.r')
    plt.xlabel('USDA')
    plt.ylabel('MODIS')
    title = 'MODIS vs. USDA cropland area [ha]'
    savefile_name = 'usda_vs_modis'
    plt.title(title)
    fig.savefig(savefile_name+'.png')
    plt.close()

# plot crop map acreage 
# regress USDA crop with MODIS
#    return acre_crop, acre_luse, fips_luse
