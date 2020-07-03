import sys, os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
#from pandas.io.json import json_normalize
import plotly.figure_factory as ff
#import plotly.express as px
#from mpl_toolkits.basemap import Basemap,maskoceans
from urllib.request import urlopen
from shapely.geometry import Point, Polygon
import json
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
from itertools import chain
from read_file import read_var_hdf
from os.path import expanduser
home = expanduser("~")

# data source: https://www.fsa.usda.gov/news-room/efoia/electronic-reading-room/frequently-requested-information/crop-acreage-data/index
# https://quickstats.nass.usda.gov/

# https://gis.stackexchange.com/questions/90055/finding-if-two-polygons-intersect-in-python
# https://automating-gis-processes.github.io/CSC18/lessons/L4/point-in-polygon.html

#  (a) Assign FIPS to each MODIS pixels: (i) read county geometry file: https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json , and (ii) check whether any pixel is located within a polygon (county) using this method: https://automating-gis-processes.github.io/CSC18/lessons/L4/point-in-polygon.html

# dealing with enclaves with Polygon: https://automating-gis-processes.github.io/site/notebooks/L1/geometric-objects.html

# for each county, extract land cover data LC_Type1 across all classes. 
# Compute surface area of cropland (and other land surface).
# load USDA estimate (convert acre in m2). Is there consistency?

# save LC_Type1 classes for all US counties, and save into file (start with 1 state)
# compare with USDA estimate. Assess whether crop and grassland estimates are coherent between USDA and MODIS datasets. Display as a scatter plot.
# For each crop pixel, load NDVI for time period, and cluster in groups. Compare group size with USDA estimates.

tickslabel = ['Evgrn Ndlf','Evgrn Brdlf','Decdu Ndlf','Decdu Brdlf','Mixd Forest','Closd Shrub','Opn Shrub','Woody Svnas','Svnas','Grasslnds','Wetlnds','Crop','Urban','Crop/Natrl','Snow/Ice','Barren','Water']

fips_list_wiki_url = "https://en.wikipedia.org/wiki/List_of_United_States_FIPS_codes_by_county"

def get_fips_list(state_name):

    data = pd.read_html(fips_list_wiki_url)
    table = data[1]
    table = np.asarray(table)

    fips = np.squeeze(table[:,0])
    name = np.squeeze(table[:,1])
    state = np.squeeze(table[:,2])
    
    if state_name!='USA':
        index = np.argwhere(state==state_name)
        fips = fips[index]
        name = name[index]
        state = state[index]

    return fips,name,state

def read_modis_landuse(state_name,varname,timestamp):

    rootpath = "/Users/xl332/Desktop/DataProjects/UsdaCrop/"
    filename = 'landuse' + '_' + timestamp + '_' + state_name + '.pkl'
    landuse_county_df = pd.read_pickle(rootpath + filename)

    return landuse_county_df

def get_modis_landuse_USA(varname,timestamp):

    state_list = ["Alabama", "Arkansas", "Arizona", "California", "Colorado", "Connecticut","Delaware", "Florida", "Georgia", "Iowa", "Idaho", "Illinois", "Indiana", "Kansas", "Kentucky", "Louisiana", "Massachusetts", "Maryland", "Maine", "Michigan", "Minnesota", "Missouri", "Mississippi", "Montana", "North Carolina", "North Dakota", "Nebraska", "New Hampshire", "New Jersey", "New Mexico", "Nevada", "New York", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Virginia", "Vermont", "Washington", "Wisconsin", "West Virginia", "Wyoming"]

    for state_name in state_list:
        print(state_name)
        get_modis_landuse(varname,timestamp,state_name)

def get_modis_landuse(varname,timestamp,state_name='USA'):

    landuse_county_df = pd.DataFrame(columns = np.concatenate((['fips'],tickslabel),axis=0))
                                    
    [fips_list,county_name_list,county_state] = get_fips_list(state_name)
    print('there are ' + str(len(fips_list)) + ' counties in ' + state_name)
    N = 0
    filename = 'landuse' + '_' + timestamp + '_' + state_name + '.pkl'    
    for fips in fips_list:
#        print(fips)
        fips = raise_exceptions(fips)
        print(str(fips[0]) + ': ' + str(N+1)  + ' out of ' + str(len(fips_list)) + ' counties')
        data_county_sum = read_modis_landuse_by_county(varname,timestamp,fips[0])
        landuse_county_df.loc[N]= np.concatenate((fips,data_county_sum),axis=0)
        N = N + 1
    landuse_county_df.to_pickle(filename)

    return landuse_county_df

def raise_exceptions(fips):

    if fips[0]==46102:
        fips=[46113]

    return fips

def read_modis_landuse_by_county(varname,timestamp,fips):

    nom_res = '500m'
    nr = get_num_pixels(nom_res)
    w = T/nr
    area_pixel_in_acre = w**2 * m2_to_acre

    [data_county_2D, lon_county_2D, lat_county_2D, extent_county] = read_modis_field_by_county(varname,timestamp,fips)

    data_county_sum = []
    for N in np.arange(1,17+1,1):
        data_county_2D_LT = np.zeros(np.shape(data_county_2D))
        data_county_2D_LT[np.where(data_county_2D==N)] = 1.0
        data_county_sum_N = np.sum(data_county_2D_LT) * area_pixel_in_acre
        data_county_sum.append(data_county_sum_N)    

    return data_county_sum

def read_modis_field_by_county(varname,timestamp,fips):

    [locator_MODIS, extent_county] = find_tiles_and_pixels(fips)

    county_pixels = locator_MODIS['pixels']
    county_H = locator_MODIS['H']  
    county_V = locator_MODIS['V']
   
    county_lon = locator_MODIS['lon']
    county_lat = locator_MODIS['lat']

###############################################            
    if isinstance(county_H,int) is True:
        ntile = 1
        H = county_H
        V = county_V
        lon = county_lon
        lat = county_lat
        pixels = county_pixels
        region = 'h'+ str('{0:02d}'.format(H)) + 'v'+ str('{0:02d}'.format(V))
        [data_tile, data_range] = read_var_hdf(varname,timestamp,region)
        data_flat = data_tile.flatten()
        data = np.nan*np.zeros(np.shape(data_flat))
        data[pixels] = data_flat[pixels]
        data_county_2D = np.reshape(data,(int(len(data)**0.5),int(len(data)**0.5)))
        lon_county_2D = np.reshape(lon,(int(len(lon)**0.5),int(len(lon)**0.5)))
        lat_county_2D = np.reshape(lat,(int(len(lat)**0.5),int(len(lat)**0.5)))
    else:
        H_min = min(county_H)
        V_min = min(county_V) 
        H_max = max(county_H)
        V_max = max(county_V)
        ntile = len(county_H)        

        tile_id = []
        for T in np.arange(0,ntile,1):
            tile_id.append(str(county_H[T])+'_'+str(county_V[T]))
        T = 0
        data_tiles_2D = []; lon_tiles_2D = []; lat_tiles_2D = []
        for V in np.arange(V_min,V_max+1,1):    
            data_tiles_H_2D = []; lon_tiles_H_2D = []; lat_tiles_H_2D = []
            for H in np.arange(H_min,H_max+1,1):
                tile_id_loc = str(H)+'_'+str(V)
                k_V = np.argwhere(V==county_V)
                k_H = np.argwhere(H==county_H)
                region = 'h'+ str('{0:02d}'.format(H)) + 'v'+ str('{0:02d}'.format(V))
                try:
                    k = np.intersect1d(k_H,k_V)[0]                    
                    lon_tile = county_lon[k]
                    lat_tile = county_lat[k]
                    pixels_tile = county_pixels[k]
                    lon_tile_2D_T = np.reshape(lon_tile,(int(len(lon_tile)**0.5),int(len(lon_tile)**0.5)))
                    lat_tile_2D_T = np.reshape(lat_tile,(int(len(lat_tile)**0.5),int(len(lat_tile)**0.5)))           
                except IndexError:
                    pixels_tile = []
                    [lon_tile_2D_T, lat_tile_2D_T] = define_SIN_axis(region,nom_res='500m')

                if len(pixels_tile)==0:
                    data_tile_2D_T = np.nan*np.zeros(np.shape(lon_tile_2D_T))
                else:
                    [data_tile, data_range] = read_var_hdf(varname,timestamp,region)
                    data_flat = data_tile.flatten()
                    d_tile = np.nan*np.zeros(np.shape(data_flat))
                    d_tile[pixels_tile] = data_flat[pixels_tile]
                    data_tile_2D_T = np.reshape(d_tile,(int(len(d_tile)**0.5),int(len(d_tile)**0.5)))

                data_tiles_H_2D.append(data_tile_2D_T)
                lon_tiles_H_2D.append(lon_tile_2D_T)
                lat_tiles_H_2D.append(lat_tile_2D_T)
            data_tiles_2D.append(data_tiles_H_2D)
            lon_tiles_2D.append(lon_tiles_H_2D)
            lat_tiles_2D.append(lat_tiles_H_2D)
 
        data_tiles_2D = np.asarray(data_tiles_2D)
        lon_tiles_2D = np.asarray(lon_tiles_2D)
        lat_tiles_2D = np.asarray(lat_tiles_2D)

        for v in np.arange(0,np.shape(data_tiles_2D)[0],1): 
            data_county_2D_H = np.squeeze(data_tiles_2D[v,0,:,:])
            lon_county_2D_H = np.squeeze(lon_tiles_2D[v,0,:,:])
            lat_county_2D_H = np.squeeze(lat_tiles_2D[v,0,:,:])
            for h in np.arange(1,np.shape(data_tiles_2D)[1],1):
                data_county_2D_H = np.concatenate((data_county_2D_H,np.squeeze(data_tiles_2D[v,h,:,:])),axis=0)
                lon_county_2D_H = np.concatenate((lon_county_2D_H,np.squeeze(lon_tiles_2D[v,h,:,:])),axis=0)
                lat_county_2D_H = np.concatenate((lat_county_2D_H,np.squeeze(lat_tiles_2D[v,h,:,:])),axis=0)
            if v==0:
                data_county_2D = data_county_2D_H
                lon_county_2D = lon_county_2D_H
                lat_county_2D = lat_county_2D_H
            else:
                data_county_2D = np.concatenate((data_county_2D_H,data_county_2D),axis=1)
                lon_county_2D = np.concatenate((lon_county_2D_H,lon_county_2D),axis=1)
                lat_county_2D = np.concatenate((lat_county_2D_H,lat_county_2D),axis=1)

    return data_county_2D, lon_county_2D, lat_county_2D, extent_county

def find_tiles_and_pixels(fips):

    if str(fips).isnumeric() is True:
        [county_shape, county_coord, county_fips] = find_county_coordinates(fips)
    else:
        [county_shape, county_coord, county_fips] = find_state_coordinates(fips)

    lon_points = []; lat_points = []; lon_points_in_holes = []; lat_points_in_holes = []
    if county_shape=='Polygon':
        if len(np.shape(county_coord))==1:
            num_enclave = np.shape(county_coord)[0]
            print('there may be ' + str(num_enclave-1) + ' enclave(s)')
            county_holes = []
            for n in np.arange(1,np.shape(county_coord)[0],1):
                county_holes = county_coord[n]
                num_points_in_holes = np.squeeze(county_holes)
                lon_points_in_holes_ = []; lat_points_in_holes_ = []
                for points in num_points_in_holes:
                    lon_points_in_holes_.append(points[0])
                    lat_points_in_holes_.append(points[1])
                lon_points_in_holes.append(lon_points_in_holes_)
                lat_points_in_holes.append(lat_points_in_holes_)
            county_coord = county_coord[0]

        num_points = np.squeeze(county_coord)
        for points in num_points:
            lon_points.append(points[0])
            lat_points.append(points[1])

    elif county_shape=='MultiPolygon':
        for n in np.arange(0,np.shape(county_coord)[0],1):            
            num_points = np.squeeze(county_coord[n])
            lon_points_p = []; lat_points_p = []
            for points in num_points:
                lon_points_p.append(points[0])
                lat_points_p.append(points[1])
            lon_points.append(lon_points_p)
            lat_points.append(lat_points_p)
            

    # locate indices in tile
    [locator_MODIS, extent_county] = find_SIN_indices(lon_points,lat_points,lon_points_in_holes,lat_points_in_holes,county_shape)

    return locator_MODIS, extent_county

def find_SIN_indices(lon_points_,lat_points_,lon_points_in_holes,lat_points_in_holes,county_shape):

    # Create a Polygon
    locator_MODIS = {}
    extent_county = {}
    H_poly_ = []; V_poly_ = []; pixels_poly_ = []; lon_poly_ = []; lat_poly_ = []
    if county_shape=='MultiPolygon':
        p = 0
        for m in np.arange(0,np.shape(lon_points_)[0],1):
            coords = []
            lon_points = lon_points_[m] 
            lat_points = lat_points_[m]
            for n in np.arange(0,len(lon_points),1):
                coords.append((lat_points[n],lon_points[n]))
            poly = Polygon(coords)
            [locator_MODIS_T, extent_county_T] = find_pixels_in_poly(lon_points,lat_points,poly)
            H_poly_.extend(locator_MODIS_T['H'])
            V_poly_.extend(locator_MODIS_T['V'])
            pixels_poly_.extend(locator_MODIS_T['pixels'])
            lon_poly_.extend(locator_MODIS_T['lon'])
            lat_poly_.extend(locator_MODIS_T['lat'])
            p = p +1
            if m==0:
                extent_county['west'] = extent_county_T['west']
                extent_county['east'] = extent_county_T['east']
                extent_county['south'] = extent_county_T['south']
                extent_county['north'] = extent_county_T['north']
            else:
                extent_county['west'] = min(extent_county['west'],extent_county_T['west'])
                extent_county['east'] = max(extent_county['east'],extent_county_T['east'])
                extent_county['south'] = min(extent_county['south'],extent_county_T['south'])
                extent_county['north'] = max(extent_county['north'],extent_county_T['north'])

    else:
        coords = []
        lon_points = lon_points_
        lat_points = lat_points_
        for n in np.arange(0,len(lon_points),1):
            coords.append((lat_points[n],lon_points[n]))

        if len(lon_points_in_holes)==0:
                poly = Polygon(coords)            
        else:
            coords_holes = []
            for m in np.arange(0,np.shape(lon_points_in_holes)[0],1):
                coords_holes_ = []
                lon_points_in_holes_ = lon_points_in_holes[m]
                lat_points_in_holes_ = lat_points_in_holes[m]
                for n in np.arange(0,len(lon_points_in_holes_),1):
                    coords_holes_.append((lat_points_in_holes_[n],lon_points_in_holes_[n]))
                coords_holes.append(coords_holes_)
            poly = Polygon(shell=coords,holes=coords_holes)

        [locator_MODIS_T, extent_county] = find_pixels_in_poly(lon_points,lat_points,poly)
        H_poly_.extend(locator_MODIS_T['H'])
        V_poly_.extend(locator_MODIS_T['V'])
        pixels_poly_.extend(locator_MODIS_T['pixels'])
        lon_poly_.extend(locator_MODIS_T['lon'])
        lat_poly_.extend(locator_MODIS_T['lat'])

    ###
    len_tile = len(H_poly_)
    tile_id_ = []; decm_id_ = []
    for T in np.arange(0,len_tile,1):
        H_T = H_poly_[T]
        V_T = V_poly_[T]
        tile_id_.append(str(H_T)+'_'+str(V_T))
        decm_id_.append(V_T+H_T*0.01)

    #### reordering must happen here
    order = np.argsort(decm_id_)
    tile_id = []; H_poly = []; V_poly = []; pixels_poly = []; lat_poly = []; lon_poly = []
    for n in np.arange(0,len(order),1):
        tile_id.append(tile_id_[order[n]])
        H_poly.append(H_poly_[order[n]])
        V_poly.append(V_poly_[order[n]])
        pixels_poly.append(pixels_poly_[order[n]])
        lat_poly.append(lat_poly_[order[n]])
        lon_poly.append(lon_poly_[order[n]])            
    ###

    tile_unique_id = []
    T = 0
    H_in_ = []; V_in_ = []; lon_in_ = []; lat_in_ = []; pixels_in_ = []
    print(tile_id)
    for tile in tile_id:
        if tile_unique_id.count(tile)==0:
            tile_unique_id.append(tile)
            H_in_.append(H_poly[T])
            V_in_.append(V_poly[T])                
            lon_in_.append(lon_poly[T])
            lat_in_.append(lat_poly[T])
            pixels_in_.append(pixels_poly[T])                 
            tile_num = tile_id.index(tile)
        else:
            tile_num = tile_unique_id.index(tile)
            pixels_in_[tile_num].extend(pixels_poly[T])
        T = T + 1

    locator_MODIS['H'] = H_in_
    locator_MODIS['V'] = V_in_
    locator_MODIS['lon'] = lon_in_
    locator_MODIS['lat'] = lat_in_
    locator_MODIS['pixels'] = pixels_in_
                
    return locator_MODIS, extent_county


def find_pixels_in_poly(lon_points,lat_points,poly):

    # locate tile     
    lon_min = min(lon_points)
    lon_max = max(lon_points)
    lat_min = min(lat_points)
    lat_max = max(lat_points)

    [H1,V1] = find_SIN_tiles(lat_max,lon_min)
    [H2,V2] = find_SIN_tiles(lat_max,lon_max)
    [H3,V3] = find_SIN_tiles(lat_min,lon_min)
    [H4,V4] = find_SIN_tiles(lat_min,lon_max)

    H_list = [H1,H2,H3,H4]
    V_list = [V1,V2,V3,V4]
    H_min = min(H_list)
    H_max = max(H_list) 
    V_min = min(V_list)
    V_max = max(V_list) 

    tiles = []
    for H in np.arange(H_min,H_max+1,1):
        for V in np.arange(V_min,V_max+1,1):
            tiles.append([H,V])

    locator_MODIS = {}
    n = 0

    MODIS_H=[]; MODIS_V=[]; MODIS_pixels=[]; MODIS_lon=[]; MODIS_lat=[]
    for tile in tiles:
        H = tile[0]
        V = tile[1]

        # NB: modify code to account for counties straddling tiles 
        # 0. define lat-lon grid within tile
        region = 'h'+ str('{0:02d}'.format(H)) + 'v'+ str('{0:02d}'.format(V))
        [lon_2d, lat_2d] = define_SIN_axis(region,nom_res='500m')
        nlon = np.shape(lon_2d)[1]
        nlat = np.shape(lat_2d)[0]

        lon_flat = lon_2d.flatten()
        lat_flat = lat_2d.flatten()

        # here detect dummy tiles by assessing whether each tile tuple actually intersect with the polygon
        # this can be done by creating another polygon from the corner coordinates of each tile
        # assess whether an intersection exists: if intersection is true, proceed as below; 
        # if not true, extract lat, lon but set pixels to some dummy value.        
        tile_north_west = (lat_2d.max(),lon_2d.min())
        tile_north_east = (lat_2d.max(),lon_2d.max()) 
        tile_south_east = (lat_2d.min(),lon_2d.max()) 
        tile_south_west = (lat_2d.min(),lon_2d.min()) 
        coords_tile = [tile_north_west,tile_north_east,tile_south_east,tile_south_west]
        Poly_tile = Polygon(coords_tile)

        if Poly_tile.intersects(poly) is True:

            # Zoom over maximum rectangle 
            i_m = np.argwhere((lon_flat>lon_min))
            i_M = np.argwhere((lon_flat<lon_max))
            j_m = np.argwhere((lat_flat>lat_min))
            j_M = np.argwhere((lat_flat<lat_max))
            is_inter = np.intersect1d(i_m,i_M)    
            js_inter = np.intersect1d(j_m,j_M)
            ks_inter = np.intersect1d(is_inter,js_inter)

            # for every pixel in pixel list, assess whether it is in the county polygon
            pixels_in_county = []
            for k in ks_inter:
                lat_k = lat_flat[k]
                lon_k = lon_flat[k]
                coord_k = Point(lat_k,lon_k)
                iswithin = coord_k.within(poly)
                if iswithin is True:
                    pixels_in_county.append(k)
                    
        else:

            pixels_in_county = []

        MODIS_H.append(H)
        MODIS_V.append(V)
        MODIS_pixels.append(pixels_in_county)
        MODIS_lon.append(lon_flat)
        MODIS_lat.append(lat_flat)

    locator_MODIS['H'] = MODIS_H
    locator_MODIS['V'] = MODIS_V
    locator_MODIS['pixels'] = MODIS_pixels
    locator_MODIS['lon'] = MODIS_lon
    locator_MODIS['lat'] = MODIS_lat

    extent_county = {}
    extent_county['west'] = lon_min
    extent_county['east'] = lon_max 
    extent_county['south'] = lat_min
    extent_county['north'] = lat_max        

    return locator_MODIS, extent_county


def find_state_coordinates(state_abbrv):

    df = pd.read_json('https://raw.githubusercontent.com/shawnbot/topogram/master/data/us-states.geojson')
    
    state_str = 'XX'
    index = 0
    while state_abbrv!=state_str:
        state_infos = df['features'][index]
        state_str = state_infos['properties']['postal']
        index = index + 1
    index = index - 1

    state_infos = df['features'][index]
    state_shape = state_infos['geometry']['type']
    state_coord = state_infos['geometry']['coordinates']
    state_strid = state_infos['properties']['code_local']

    state_id = int(''.join([n for n in state_strid if n.isdigit()]))

    return state_shape, state_coord, state_id

def find_county_coordinates(fips):

    # Using fips, find county coordinates
    df = pd.read_json('https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json')

    county_fips = 0
    index = 0
    while fips!=int(county_fips):
        county_infos = df['features'][index]
        county_fips = county_infos['id']
        index = index + 1
    index = index - 1

    county_infos = df['features'][index]
#    print(county_infos)
    county_shape = county_infos['geometry']['type']
    county_coord = county_infos['geometry']['coordinates']
#    county_name = county_infos['properties']['NAME']
#    county_area = county_infos['properties']['CENSUSAREA']
    county_fips = county_infos['id']
    
    return county_shape, county_coord, county_fips


def find_SIN_tiles(lat,lon,nom_res='500m'):

    nr = get_num_pixels(nom_res)
    w = T / nr

    lon_rad = lon * np.pi/180
    lat_rad = lat * np.pi/180
    x = R * lon_rad * np.cos(lat_rad)
    y = R * lat_rad 

    H = int(np.floor((x - xmin)/T))
    V = int(np.floor((ymax - y)/T))

    return H, V

def define_SIN_axis(region,nom_res='500m'):

    # the actual size of a "X-m" MODIS sinusoidal grid cell [m]
    nr = get_num_pixels(nom_res)
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

def get_num_pixels(nom_res):

    if nom_res=='500m':
        nr = 2400
    elif nom_res=='250m':
        nr = 4800
    elif nom_res=='1000m':
        nr = 1200

    return nr


# the radius of the idealized sphere representing the Earth [m]
R = 6371007.181
# the height and width of each MODIS tile in the projection plane [m]
T = 1111950.0
# the western limit of the projection plane [m]
xmin = -20015109.0
# the northern limit of the projection plane [m]
ymax = 10007555.0
# the actual size of a "X-m" MODIS sinusoidal grid cell [m]    
# conversion factor from m2 to acre
#m2_to_acre = 0.00024710538146717
m2_to_acre = 0.0001
################# PLOTTING #################
def plot_data_county_2D(varname,timestamp,fips):

    [data_county_2D, lon_county_2D, lat_county_2D, extent_county] = read_modis_field_by_county(varname,timestamp,fips)

#    PLOT
#    print(np.shape(data_county),np.shape(lon),np.shape(lat))
    rivers_50m = cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', '50m')
    delta_deg = 0.02
    colors = ['#05450a', '#086a10', '#54a708', '#78d203', '#009900', '#c6b044', '#dcd159', '#dade48', '#fbff13', '#b6ff05', '#27ff87', '#c24f44', '#a5a5a5', '#ff6d4c', '#69fff8', '#f9ffa4', '#1c0dff']

    cm = mcolors.LinearSegmentedColormap.from_list('LandUse',colors,N=17)

    fig = plt.figure()
    extent = [extent_county['west']-delta_deg, extent_county['east']+delta_deg, extent_county['south']-delta_deg, extent_county['north']+delta_deg]    
    lon_center = np.mean([extent_county['west'], extent_county['east']])
    lat_center = np.mean([extent_county['south'], extent_county['north']])
    ax = plt.axes(projection=ccrs.PlateCarree())
#    ax = plt.axes(projection=ccrs.AlbersEqualArea())
#    rotated_pole = ccrs.RotatedPole(pole_longitude=lon_center, pole_latitude=lat_center)
    psm = ax.pcolormesh(lon_county_2D, lat_county_2D, data_county_2D,cmap=cm, vmin=1, vmax=17) #, transform=rotated_pole)
    ax.gridlines()
    ax.set_extent(extent)    
    cbar = plt.colorbar(psm, ax=ax, ticks=np.arange(1,17+1,1), fraction=0.046, pad=0.04)
    cbar.ax.set_yticklabels(tickslabel)
    cbar.ax.tick_params(labelsize=6)

#    ax.coastlines(resolution='50m') 
#    ax.coastlines()
#    ax.add_feature(cartopy.feature.RIVERS)
    ax.add_feature(rivers_50m, facecolor='None', edgecolor='b')
#    ax.add_feature(cartopy.feature.COASTLINE)
#    ax.add_feature(cartopy.feature.OCEAN)
#    ax.add_feature(cartopy.feature.LAND, edgecolor='black')
    save_file = 'tmp.png'
    fig.savefig(save_file)
    plt.close()

    # plot on local pixels on map: lon_2d, lat_2d, var [opt]
    # plot map of counties on pixelized maps [opt]
    # perform budget of land use (cropland+grassland?) and plot on county map
    # compare crop acreage by county from USDA and MODIS 
#    ax = plt.axes(projection=ccrs.Orthographic(lon_center,lat_center))
    # Start statistics with Vegetation index
