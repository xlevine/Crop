import numpy as np
from os import listdir
from os.path import isfile, join
from pyhdf.SD import SD, SDC
from dic_modis import find_dataset_varname, data_folder, grid_coefficient

## root path must be changed to the user's tree data directory
rootpath = "/Users/xl332/Desktop/DataProjects/LandCloud/Data/"

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
        
    [datafolder, time_format, nom_res]= data_folder(dataname,timestamp)
    spatial_stamp = grid_coefficient(dataname,region)

    year = int(str(timestamp)[:4])
    folders_path = rootpath + datafolder + "/"

    hdf_list = [f for f in listdir(folders_path) if isfile(join(folders_path, f))]
    hdf_prec = datafolder + ".A" + str(year) + time_format + spatial_stamp
    file_names = [i for i in hdf_list if i.startswith(hdf_prec)]    
    file_name = file_names[0]
    filename = folders_path + file_name

    return filename

def read_var_hdf(varname,timestamp,region='global'):

    dataname = find_dataset_varname(varname,timestamp)
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
 
    if '_SIN' in dataname:
        H = int(region[1:3])
        V = int(region[4:]) 
        data = data.T
        if H > 17:
            data = np.flip(data,axis=0)
    else:
        data = np.flip(data,axis=0)

    return data, data_range
