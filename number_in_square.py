# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 10:31:18 2015

@author: zhaobin
"""
'plot number of observation in 1 minute square'

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np
import sys
sys.path.append('../modules')
from conversions_old import dm2dd
from turtleModule import draw_basemap,whichArea
########################################################################
one_minute=1.0/60   #1.0/60 is 1 minute
data=pd.read_csv('binned_td_hoey.csv')
lat,lon=dm2dd(data['LATITUDE'],data['LONGITUDE']) #convert ddmm.m to dd.ddd
lat_i=np.arange(min(lat),max(lat),one_minute)   
lon_i=np.arange(min(lon),max(lon),one_minute)    
number=[]
for i in range(len(lon_i)):
    j=[0]*len(lat_i)
    number.append(j)
lat_n=[]
lon_m=[]
for i in range(len(lat)):
    n=whichArea(lat[i],lat_i)
    m=whichArea(lon[i],lon_i)
    number[m][n]+=1                                      # calculate number in 1 minute square
lonsize=[min(lon),max(lon)]
latsize=[min(lat),max(lat)]
lon_is = np.linspace(lonsize[0],lonsize[1],1000)
lat_is = np.linspace(latsize[0],latsize[1],1000)     
LON,LAT,NUMBER=[],[],[]
for j in range(len(lon_i)):
    for i in range(len(lat_i)):
        LON.append(lon_i[j])
        LAT.append(lat_i[i])
        NUMBER.append(number[j][i])                    #use for griddata
number_i = griddata(np.array(LON),np.array(LAT),np.array(NUMBER),lon_is,lat_is,interp='linear')
depth_i = griddata(np.array(lon),np.array(lat),np.array(data['MEAN_DEPTH']),lon_is,lat_is,interp='linear')

fig=plt.figure(figsize=(12,10))
ax=fig.add_subplot(111)
draw_basemap(fig,ax,lonsize,latsize)
#CS1=plt.contour(lon_is, lat_is,depth_i,np.arange(100,102,100),colors='r',linewidths=4,linestyles=':')   #plot 100m depth
CS = plt.contourf(lon_is, lat_is, number_i, np.arange(1,101,1), cmap=plt.cm.rainbow,vmax=101, vmin=1)
cbar=plt.colorbar(CS,ticks=[1,20,40,60,80,100])
cbar.ax.tick_params(labelsize=20)
#x=(min(lon)+176*one_minute,min(lon)+430*one_minute,min(lon)+236*one_minute,min(lon)+288*one_minute)
#y=(min(lat)+254*one_minute,min(lat)+300*one_minute,min(lat)+252*one_minute,min(lat)+406*one_minute)
#plt.scatter(x,y,c='y',s=50)
plt.title('Number of observation in 1 minute square(2006~2015)',fontsize=20)
plt.show()
