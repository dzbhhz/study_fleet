# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 10:43:32 2015

@author: jmanning
"""
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
import matplotlib as mpl
import math
import numpy as np
from conversions_old import dm2dd
from turtleModule import whichArea,dist
def get_roms(lat,lon,depth,time):
    url_1='http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/hidden/2006_da/his'
    nc_1=netCDF4.Dataset(url_1)
    time_1=nc_1.variables['ocean_time'][:]
    lats=nc_1.variables['lat_rho'][:]
    lons=nc_1.variables['lon_rho'][:]
    s_rho=nc_1.variables['s_rho'][:]
    h=nc_1.variables['h'][:]
    temp_1=nc_1.variables['temp']
    url_2='http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/his_Best/ESPRESSO_Real-Time_v2_History_Best_Available_best.ncd'
    nc_2=netCDF4.Dataset(url_2)
    time_2=nc_2.variables['time'][:]
    temp_2=nc_2.variables['temp']
    TEMP=[]
    for i in depth.index:
        Distance=dist(lat[i],lon[i],lats,lons)
        near_index=np.argmin(Distance)
        index_one,index_two=near_index/len(Distance[0]),near_index%len(Distance[0])
        H=h[index_one,index_two]
        layer_depth=-s_rho*H
        index_D=np.argmin(np.abs(depth[i]-layer_depth))
        if (time[i]-datetime(2013,05,18)).total_seconds()/3600>25: #25 hour is the first time ni 2013_da
            index_T=np.argmin(np.abs((time[i]-datetime(2013,05,18)).total_seconds()/3600-time_2))
            tt=temp_2[index_T,index_D,index_one,index_two]
            TEMP.append(tt)
        else:
            index_T=np.argmin(np.abs((time[i]-datetime(2006,01,01)).total_seconds()-time_1))
            tt=temp_1[index_T,index_D,index_one,index_two]
            TEMP.append(tt)
        print i
    return TEMP
#####################################################################
one_minute=1.0/60   #1.0/60 is 1 minute
data=pd.read_csv('binned_td_FVCOM.csv')
lat,lon=dm2dd(data['LATITUDE'],data['LONGITUDE']) #convert ddmm.m to dd.ddd
for i in data.index:
    data['datet'][i]=datetime.strptime(data['datet'][i],'%Y-%m-%d %H:%M:%S')
temp_roms=get_roms(lat,lon,data['MEAN_DEPTH'],data['datet'])
for i in range(len(temp_roms)):
    if temp_roms[i]<100:
        temp_roms[i]=temp_roms[i]
    else:
        temp_roms[i]=-10
data['temp_roms']=pd.Series(temp_roms)
data.to_csv('binned_td_FVCOM_ROMS.csv')
