# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:41:15 2015
plot depth of site 1 and site 2
@author: zhobin
"""
import numpy as np 
import pandas as pd
import matplotlib
import netCDF4 
import datetime
import pytz
import matplotlib.pyplot as plt
from MODULES import dm2dd,DIST,FIG_PLOT_emolt,whichArea,draw_basemap,stick_plot,get_emolt_data,getobs_tempsalt
from pylab import *
import matplotlib.tri as Tri
from matplotlib.mlab import griddata
one_minute=1.0/60
sites=[[185,281],[193,273]]   #different sites to plot
data=pd.read_csv('binned_td_FVCOM.csv')
lat,lon=dm2dd(data['LATITUDE'],data['LONGITUDE']) #convert ddmm.m to dd.ddd
lat_i=np.arange(min(lat),max(lat),one_minute)   
lon_i=np.arange(min(lon),max(lon),one_minute)
alat=lat_i[sites[0][1]]
alon=lon_i[sites[0][0]]
lonsize=[alon-20*one_minute,alon+20*one_minute]
latsize=[alat-20*one_minute,alat+20*one_minute]
url='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'
nc=netCDF4.Dataset(url)
lons=nc.variables['lon'][:]
lats=nc.variables['lat'][:] 
h = nc.variables['h'][:] 
tri = Tri.Triangulation(lons,lats) 
indx=argwhere((lons >= lonsize[0]) & (lons <= lonsize[1]) & (lats >= latsize[0]) & (lats <= latsize[1]))  # find  points in bounding box
#plot picture
fig=plt.figure(figsize=(12,10))
ax=fig.add_subplot(111)
draw_basemap(fig, ax, lonsize, latsize,interval_lon=10*one_minute, interval_lat=10*one_minute)
levels=arange(int(-max(h[indx])),0.1,0.5)
tricontourf(tri,-h,levels=levels,shading='faceted',cmap=plt.cm.rainbow)
gca().patch.set_facecolor('0.5')
cbar=colorbar()
cbar.set_label('Water Depth (m)')
ax.scatter(lon_i[sites[0][0]],lat_i[sites[0][1]],s=200,c='red',marker='*')
ax.scatter(lon_i[sites[1][0]],lat_i[sites[1][1]],s=200,c='red',marker='*')
plt.annotate('site_1', xy=(lon_i[sites[0][0]], lat_i[sites[0][1]]),xytext=(lon_i[sites[0][0]],lat_i[sites[0][1]]-3*one_minute),arrowprops=dict(arrowstyle='->'),fontsize=12)
plt.annotate('site_2', xy=(lon_i[sites[1][0]], lat_i[sites[1][1]]),xytext=(lon_i[sites[1][0]],lat_i[sites[1][1]]-3*one_minute),arrowprops=dict(arrowstyle='->'),fontsize=12)
plt.title('Depth of site 1 and site 2',fontsize=20)
plt.show()
