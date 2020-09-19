def data_folder(dataset_name,timestamp):

    if dataset_name=="BurnedArea_SIN":
        # 500 m SIN monthly
        datafolder = "MCD64A1"    
        year = int(str(timestamp)[:4])
        month_num = int(str(timestamp)[-2:])
        first_day_in_month = first_day_month(year,month_num)
        time_format = '{0:03d}'.format(first_day_in_month)
        nom_res = '500m' 
    elif dataset_name=="LandCover_SIN":
        # 500 m SIN yearly
        datafolder = "MCD12Q1"    
        time_format = "001"
        nom_res = '500m' 
    elif dataset_name=="LandDyn_SIN":
        # 500 m SIN yearly
        datafolder = "MCD12Q2"    
        time_format = "001"
        nom_res = '500m' 
    elif dataset_name=="VegIndex_SIN":
        # 500 m SIN yearly
        datafolder = "MOD13A1"    
        day_out_16 = day_out_16(year,month_num)
        time_format = '{0:03d}'.format(day_out_16)
        nom_res = '500m' 

    return datafolder, time_format, nom_res

def grid_coefficient(dataset_name,region='global'):

    if '_SIN' in dataset_name:
        tile_number = region
        spatial_stamp = "." + tile_number + ".006"
    else:
        spatial_stamp = ".006"

    return spatial_stamp

def find_dataset_varname(varname,timestamp):

    # using timestamp may be necessary to different between datasets with same variables but output at different frequencies
    # add more datasets and define variable list in DATASETNAME_GRID_list  
    if varname in BurnedArea_SIN_list:
        dataname = 'BurnedArea_SIN'
    elif varname in LandDyn_SIN_list:
        dataname = 'LandDyn_SIN'
    elif varname in LandCover_SIN_list:
        dataname = 'LandCover_SIN'

    return dataname

############################# VARIABLES LIST FOR EACH DATASET (can be found in MODIS portal description of each dataset)

BurnedArea_SIN_list = ['Burn Date', 'First Day']
LandDyn_SIN_list = ['Dormancy', 'MidGreendown', 'NumCycles', 'Greenup', 'Senescence', 'EVI_Amplitude', 'EVI_Area', 'Maturity', 'MidGreenup', 'EVI_Minimum', 'Peak']
LandCover_SIN_list = ['QC', 'LC_Type4', 'LC_Type5', 'LC_Prop1_Assessment', 'LC_Type3', 'LC_Type1', 'LC_Prop2_Assessment', 'LW', 'LC_Prop3_Assessment', 'LC_Prop3', 'LC_Prop2', 'LC_Prop1', 'LC_Type2']
