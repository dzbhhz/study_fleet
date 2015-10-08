# -*- coding: utf-8 -*-
"""
Created on Thu May  7 09:59:45 2015
1)write a file named binned_td_FVCOM.csv which has FVCOM temperature
2)write a file named binned_td_FVCOM_ROMS.csv which has ROMS temperature
@author: zhaobin
"""
import numpy as np
import pandas as pd
import netCDF4
from datetime import datetime,timedelta
from conversions_old import dm2dd
from turtleModule import draw_basemap,whichArea,dist,str2ndlist
def getFVcom(latpt,lonpt,time_roms):
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
    nc = netCDF4.Dataset(url)
    lats = nc.variables['lat'][:]
    lons = nc.variables['lon'][:]
    Depth = nc.variables['h'][:]
    siglay=nc.variables['siglay'][:]
    time=nc.variables['time'][200000:316008]#time between 2006-2013
    modtemp=nc.variables['temp'] 
    modtime=[]
    for i in range(len(time)):
        t=timedelta(days=float(time[i])).total_seconds()
        modtime.append(t)                        #change time days to seconds
    depth=(-Depth*siglay).transpose()                #each layer`s depth
    Temp=[]
    for i in latpt.index:
        distance=dist(lonpt[i],latpt[i],lons,lats)
        node=np.argmin(distance)                 #get nearest node
        t_diff=(time_roms[i]-datetime(1858,11,17)).total_seconds()     #1858,11,17 is FVCOM`s start time
        TIME=np.argmin(abs(t_diff-np.array(modtime)))
        T=modtemp[200000+TIME][44][node]
        Temp.append(T)
    nc.close()
    return Temp
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
            if type(tt) not float:
                tt=-10
            TEMP.append(tt)
        else:
            index_T=np.argmin(np.abs((time[i]-datetime(2006,01,01)).total_seconds()-time_1))
            tt=temp_1[index_T,index_D,index_one,index_two]
            if type(tt) not float:
                tt=-10
            TEMP.append(tt)
        print i
    return TEMP
###########################################################
data=pd.read_csv('binned_td_hoey.csv')
for i in range(len(data)):
    data['datet'][i]=datetime.strptime(data['datet'][i],'%Y-%m-%d %H:%M:%S')
    if data['LATITUDE'][i]%100/60==1:
        data['LATITUDE'][i]=(data['LATITUDE'][i]/100+1)*100
    if data['LONGITUDE'][i]%100/60==1:
        data['LONGITUDE'][i]=(data['LONGITUDE'][i]/100+1)*100    #for example:change 4260 to 4300 
lat,lon=dm2dd(data['LATITUDE'],data['LONGITUDE']) #convert ddmm.m to dd.ddd
lat,lon=pd.Series(lat),pd.Series(lon)
index=[]
for i in range(len(data)):
    a=np.where(data['LATITUDE']==data['LATITUDE'][i])[0]
    b=np.where(data['LONGITUDE']==data['LONGITUDE'][i])[0]
    index.append(set(a)&set(b))
unique_index = [list(x) for x in set(tuple(x) for x in index)]
temp=[]

for i in range(len(unique_index)):
    t=getFVcom(lat[unique_index[i]],lon[unique_index[i]],data['datet'][unique_index[i]])
    temp.append(t)
temp=pd.Series(temp)

for i in range(len(data)):
    for j in range(len(unique_index)):
        for q in range(len(unique_index[j])):
            if i==unique_index[j][q]:
                Temp.append(temp[j][q])

data['temp_fvcom']=pd.Series(Temp)
data.to_csv('binned_td_FVCOM.csv')

data=pd.read_csv('binned_td_FVCOM.csv')
lat,lon=dm2dd(data['LATITUDE'],data['LONGITUDE']) #convert ddmm.m to dd.ddd
for i in data.index:
    data['datet'][i]=datetime.strptime(data['datet'][i],'%Y-%m-%d %H:%M:%S')
temp_roms=get_roms(lat,lon,data['MEAN_DEPTH'],data['datet'])
for i in range(len(temp_roms)):   #get rid of some bad data
    if temp_roms[i]<100:
        temp_roms[i]=temp_roms[i]
    else:
        temp_roms[i]=-10  
data['temp_roms']=pd.Series(temp_roms)
data.to_csv('binned_td_FVCOM_ROMS.csv')        
