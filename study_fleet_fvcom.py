# -*- coding: utf-8 -*-
"""
Created on Thu May  7 09:59:45 2015

@author: zhaobin
"""
'write a file named binned_td_fvcom.csv which has FVCOM temperature'
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
    print i 
#DATA=pd.read_csv('study_temp.csv',skiprows=0,names=['index','temp'])
temp=pd.Series(temp)
Temp=[]
for i in range(len(data)):
    for j in range(len(unique_index)):
        for q in range(len(unique_index[j])):
            if i==unique_index[j][q]:
                Temp.append(temp[j][q])
    print i
data['temp_fvcom']=pd.Series(Temp)
data.to_csv('binned_td_fvcom.csv')
          
