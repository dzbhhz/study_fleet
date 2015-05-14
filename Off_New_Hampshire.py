# -*- coding: utf-8 -*-
"""
Created on Fri May  1 11:10:45 2015

@author: zhobin
"""
'plot temperature off New Hampshire and some details'
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation
from datetime import datetime,timedelta
import netCDF4
import numpy as np
from conversions_old import dm2dd
from turtleModule import draw_basemap,dist,whichArea
#######################################################################################
one_minute=1.0/60
data=pd.read_csv('binned_td_fvcom.csv')
temp_fvcom=pd.Series(data['temp_fvcom'])
for i in range(len(data)):
    data['datet'][i]=datetime.strptime(data['datet'][i], "%Y-%m-%d %H:%M:%S")
lat,lon=dm2dd(data['LATITUDE'],data['LONGITUDE']) #convert ddmm.m to dd.ddd
lon=pd.Series(lon)
lat=pd.Series(lat)
lat_i=np.arange(min(lat),max(lat),one_minute)   
lon_i=np.arange(min(lon),max(lon),one_minute)    #1.0/60 is 1 minute
indx=[]
for i in range(len(lat)):
    n=whichArea(lat[i],lat_i)
    m=whichArea(lon[i],lon_i)
    if n==370 and m==245:
        indx.append(i)

time=[2009,2011,2012]
month=[[11,5,17],[6,4,116],[6,18,30]]
Temp,INDX_month=[[],[],[]],[[],[],[]]
for i in range(len(time)):
    for j in range(len(indx)):
        if data['datet'][indx[j]].year==time[i]:
            if data['datet'][indx[j]].month==month[i][0]:
                if month[i][1]<=data['datet'][indx[j]].day<=month[i][-1]:
                    INDX_month[i].append(indx[j])
                    Temp[i].append(j)

fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot(data['datet'][indx].values,data['MEAN_TEMP'][indx].values,label='obs')
plt.plot(data['datet'][indx].values,data['temp_fvcom'][indx].values,label='mod')
plt.legend(loc='best')
plt.ylabel('Temperature (degC)',fontsize=15)
plt.xlabel('Time',fontsize=15)
plt.title('Temparature vs time off New Hampshire',fontsize=20)

fig=plt.figure()
ax1=fig.add_subplot(311)
ax1.plot(data['datet'][INDX_month[0]].values,data['MEAN_TEMP'][INDX_month[0]].values,label='obs',linewidth=5,c='black')
ax1.plot(data['datet'][INDX_month[0]].values,temp_fvcom[Temp[0]],linestyle=':',label='mod',linewidth=5,c='black')
ax1.grid(True)
ax2=fig.add_subplot(312)
ax2.plot(data['datet'][INDX_month[1]].values,data['MEAN_TEMP'][INDX_month[1]].values,linewidth=5,c='black')
ax2.plot(data['datet'][INDX_month[1]].values,temp_fvcom[Temp[1]],linestyle=':',linewidth=5,c='black')
ax2.grid(True)
ax3=fig.add_subplot(313)
ax3.plot(data['datet'][INDX_month[2]].values,data['MEAN_TEMP'][INDX_month[2]].values,linewidth=5,c='black')
ax3.plot(data['datet'][INDX_month[2]].values,temp_fvcom[Temp[2]],linestyle=':',linewidth=5,c='black')
ax3.grid(True)
ax1.legend(loc='best')
ax1.set_title('Temperature vs time off New Hampshire',fontsize=20)
ax2.set_ylabel('Temperature(degC)',fontsize=25)
ax3.set_xlabel('Location:'+str(round(lat[indx[0]],2))+','+str(round(lon[indx[0]],2)),fontsize=25)
plt.show()
