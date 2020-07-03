import numpy as np
from os import listdir
from os.path import isfile, join
from pyhdf.SD import SD, SDC
import pprint

#datasets = ["VegIndex","LandCover"]
datasets = ["BurnedArea_SIN","LandCover_SIN","LandDyn_SIN"]
SIN_list = ["BurnedArea_SIN","LandCover_SIN","LandDyn_SIN"]
CMG_list = ["VegIndex","LandCover"]
# SDS_NAME = "Majority_Land_Cover_Type_1"
# BurnedArea_SIN = 'Burn Date', 'First Day'
# LandDyn_SIN = 'Dormancy', 'MidGreendown', 'NumCycles', 'Greenup', 'Senescence', 'EVI_Amplitude', 'EVI_Area', 'Maturity', 'MidGreenup', 'EVI_Minimum', 'Peak'
# LandCover_SIN = 'QC', 'LC_Type4', 'LC_Type5', 'LC_Prop1_Assessment', 'LC_Type3', 'LC_Type1', 'LC_Prop2_Assessment', 'LW', 'LC_Prop3_Assessment', 'LC_Prop3', 'LC_Prop2', 'LC_Prop1', 'LC_Type2'

def data_folder(dataset_name,timestamp):

#    if dataset_name=='VegIndex':
#        # 0.05deg CMG monthly        
#        datafolder = "MOD13C2"
#        year = int(str(timestamp)[:4])
#        month_num = int(str(timestamp)[-2:])
#        first_day_in_month = first_day_month(year,month_num)
#        time_format = '{0:03d}'.format(first_day_in_month)
#    elif dataset_name=='LandCover':
#        # 0.05deg CMG yearly
#        datafolder = "MCD12C1"    
#        time_format = "001"

    if dataset_name=="BurnedArea_SIN":
        # 500 m SIN monthly
        datafolder = "MCD64A1"    
        year = int(str(timestamp)[:4])
        month_num = int(str(timestamp)[-2:])
        first_day_in_month = first_day_month(year,month_num)
        time_format = '{0:03d}'.format(first_day_in_month)
    elif dataset_name=="LandCover_SIN":
        # 500 m SIN yearly
        datafolder = "MCD12Q1"    
        time_format = "001"
    elif dataset_name=="LandDyn_SIN":
        # 500 m SIN yearly
        datafolder = "MCD12Q2"    
        time_format = "001"

    elif dataset_name=="VegIndex_SIN":
        # 500 m SIN yearly
        datafolder = "MOD13A1"    
        day_out_16 = day_out_16(year,month_num)
        time_format = '{0:03d}'.format(day_out_16)

    elif dataset_name=="CloudProp_M":
        # 1.0deg X monthly
        datafolder = "CLDPROP_M3_MODIS_Aqua"    
    elif dataset_name=="CloudProp_D":
        # 1.0deg X daily
        datafolder = "CLDPROP_D3_MODIS_Aqua"    

    return datafolder, time_format

def grid_coefficient(dataset_name,region):

    if dataset_name=='LandCover' or dataset_name== 'VegIndex':
        spatial_stamp = ".006"
    else:
        if region=='global':
            region = 'h11v03'
        tile_number = region
        spatial_stamp = "." + tile_number + ".006"

    return spatial_stamp

def first_day_month(year,month_num):

    if month_num == 1:
        first_day_month = 1
    else:
        d = 1
        for m in np.arange(1,month_num,1):
            days_month = days_in_month(year,m)
            d = d + days_month
        first_day_month = d

    return first_day_month

def days_in_month(year,month):

    if month == 9 or month == 4 or month == 6 or month == 11:
        days_month = 30
    elif month == 1 or month == 3 or month == 5 or month== 7 or month == 8 or month == 10 or month== 12:
        days_month = 31
    elif month == 2 and is_leap_year(year) == True:
        days_month = 29
    elif month == 2 and is_leap_year(year) == False:
        days_month = 28
    else:
        print('ERROR')

    return days_month

def is_leap_year(year):

    return (year % 4 == 0) and (year % 100 != 0) or (year % 400 == 0)

def read_filename_hdf(dataname,timestamp,region='global'):
        
    [datafolder, time_format]= data_folder(dataname,timestamp)
    spatial_stamp = grid_coefficient(dataname,region)

    year = int(str(timestamp)[:4])
    folders_path =  "/Users/xl332/Desktop/DataProjects/LandCloud/Data/" + datafolder + "/"

    hdf_list = [f for f in listdir(folders_path) if isfile(join(folders_path, f))]
    hdf_prec = datafolder + ".A" + str(year) + time_format + spatial_stamp
    file_names = [i for i in hdf_list if i.startswith(hdf_prec)]    
    file_name = file_names[0]
    filename = folders_path + file_name

    return filename

def get_varlist_hdf(dataname,timestamp,region='global'):

    if len(str(timestamp))==4:
        timestamp = timestamp + '01'

    filename = read_filename_hdf(dataname,timestamp,region)
    file = SD(filename,SDC.READ)
    datasets_dic = file.datasets()
    sds_list = []
    for idx,sds in enumerate(datasets_dic.keys()):
        sds_list.append(sds)
#        print idx,sds

    return sds_list

def find_dataset_varname(varname,timestamp,region='global'):
    
#    for dataset in datasets:
#        sds_list = get_varlist_hdf(dataset,timestamp,region)            
#        if any(varname in s for s in sds_list):
#            dataname = dataset
    if varname=='LC_Type1':
        dataname = 'LandCover_SIN'

    return dataname

def read_var_hdf(varname,timestamp,region='global'):

    dataname = find_dataset_varname(varname,timestamp,region)
    filename = read_filename_hdf(dataname,timestamp,region)

    file = SD(filename,SDC.READ)
 
    # select and read the sds data
    sds = file.select(varname)
    data = sds.get()

    add_offset=0; scale_factor=1.0
    for key in sds.attributes().keys():
        if key == 'add_offset':
            add_offset = sds.add_offset  
        if key == 'scale_factor':
            scale_factor = sds.scale_factor
        if key == '_FillValue':
            fill_value = sds._FillValue
        if key == 'valid_range':
            valid_range = sds.valid_range

    data = (data - add_offset) * scale_factor
    fill_value = (fill_value - add_offset) * scale_factor
    data[np.where(data==fill_value)] = np.nan
    max_range = (valid_range[1] - add_offset) * scale_factor
    min_range = (valid_range[0] - add_offset) * scale_factor
    data_range = [min_range, max_range]
 
    if dataname in SIN_list:
        H = int(region[1:3])
        V = int(region[4:]) 
        data = data.T
        if H > 17:
            data = np.flip(data,axis=0)
    elif dataname in CMG_list:
        data = np.flip(data,axis=0)

    return data, data_range

def define_axis(varname,timestamp,region='global'):

    dataname = find_dataset_varname(varname,timestamp,region)
    if dataname in SIN_list:
        [lon_2d, lat_2d] = define_SIN_axis(region)
    elif dataname in CMG_list:
        [lon_2d, lat_2d] = define_CMG_axis()

    return lon_2d, lat_2d

def define_CMG_axis():

    delta = 0.05
    lat = np.arange(-90.00+delta/2,90.00,delta)
    lon = np.arange(-180.00+delta/2,180.00,delta)
    nlon = len(lon); nlat = len(lat)

    lat_2d = np.tile(lat,(nlon,1)).T
    lon_2d = np.tile(lon,(nlat,1))

    return lon_2d, lat_2d

def define_SIN_axis(region,nom_res='500m'):

    # the radius of the idealized sphere representing the Earth [m]
    R = 6371007.181
    # the height and width of each MODIS tile in the projection plane [m]
    T = 1111950.0
    # the western limit of the projection plane [m]
    xmin = -20015109.0
    # the northern limit of the projection plane [m]
    ymax = 10007555.0
    # the actual size of a "X-m" MODIS sinusoidal grid cell [m]
    if nom_res=='500m':
        nr = 2400
    elif nom_res=='250m':
        nr = 4800
    elif nom_res=='1000m':
        nr = 1200
    w = T / nr

    H = int(region[1:3])
    V = int(region[4:])
    i = np.arange(0,nr,1)
    j = np.arange(0,nr,1)    
    x = (j + 0.5)*w + H*T + xmin
    y = ymax - (i + 0.5)*w - V*T

    lat_2d = np.zeros((nr,nr)); lon_2d = np.zeros((nr,nr))

    for indy in i:
        lat_2d[indy,:] = y[indy] / R
        lon_2d[indy,:] = np.divide(x[:],(R * np.cos(y[indy] / R)))

    lat_2d = np.asarray(lat_2d); lon_2d = np.asarray(lon_2d)

    lat_2d = lat_2d * 180.0/np.pi; lon_2d = lon_2d * 180.0/np.pi
    lat_2d = lat_2d.T
    lon_2d = lon_2d.T
    if H > 17:
        lon_2d = np.flip(lon_2d,axis=0)    

    lat_2d = np.flip(lat_2d,axis=0)

    return lon_2d, lat_2d

# The MODIS Climate Modeling Grid (CMG) consists of 7200 columns by 3600 rows. Each cell has a resolution of 0.05 (approximately 5 km). The upper-left corner of the upper-left cell is -180.00 longitude, 90.00 latitude. The lower-right corner of the lower right cell is -180.00 longitude, -90.00 latitude. For additional details about the MODIS Climate Modeling Grid, see the NASA MODIS Lands | MODIS Grids Web page. https://modis-land.gsfc.nasa.gov/MODLAND_grid.html

############################# 
