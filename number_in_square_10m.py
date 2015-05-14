# -*- coding: utf-8 -*-
"""
Created on Fri May  8 10:42:16 2015

@author: zhaobin
"""
'plot number of observation in 10 minute square'
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
from datetime import datetime,timedelta
import numpy as np
from conversions_old import dm2dd
from turtleModule import draw_basemap,whichArea
###############################################################
ten_minute=10.0/60
data=pd.read_csv('binned_td_hoey.csv')
lat,lon=dm2dd(data['LATITUDE'],data['LONGITUDE']) #convert ddmm.m to dd.ddd
for i in range(len(data)):
    data['datet'][i]=datetime.strptime(data['datet'][i],'%Y-%m-%d %H:%M:%S')
lat_i=np.arange(int(min(lat)),int(max(lat))+1,ten_minute)   
lon_i=np.arange(int(min(lon)),int(max(lon))+1,ten_minute) 
number=[]
for q in range(10):
    number.append([])
    for i in range(len(lon_i)):
        j=[0]*len(lat_i)
        number[q].append(j)
lat_n=[]
lon_m=[]
for i in range(len(lat)):
    for j in range(10):                  #10 is 10 years,from 2006 to 2015
        if data['datet'][i].year==2006+j:
            n=whichArea(lat[i],lat_i)
            m=whichArea(lon[i],lon_i)
            number[j][m][n]+=1                                      # calculate number in 10 minute square
lonsize=[min(lon),max(lon)]
latsize=[min(lat),max(lat)]

for q in range(len(number)):
    fig=plt.figure(figsize=(12,10))
    ax=fig.add_subplot(111)
    draw_basemap(fig,ax,lonsize,latsize)
    for i in range(len(number[q])):
        for j in range(len(number[q][i])):
            if 0<number[q][i][j]<50:
                lat_10m=int(min(lat))+ten_minute*j
                lon_10m=int(min(lon))+ten_minute*i
                ax.add_patch(patches.Rectangle((lon_10m, lat_10m),ten_minute,ten_minute,color='blue'))
                p1=Rectangle((0, 0), 1, 1, fc="blue")
            if 50<=number[q][i][j]<100:
                lat_10m=int(min(lat))+ten_minute*j
                lon_10m=int(min(lon))+ten_minute*i
                ax.add_patch(patches.Rectangle((lon_10m, lat_10m),ten_minute,ten_minute,color='yellow'))
                p2=Rectangle((0, 0), 1, 1, fc="yellow")
            if 100<=number[q][i][j]<150:
                lat_10m=int(min(lat))+ten_minute*j
                lon_10m=int(min(lon))+ten_minute*i
                ax.add_patch(patches.Rectangle((lon_10m, lat_10m),ten_minute,ten_minute,color='orange'))
                p3=Rectangle((0, 0), 1, 1, fc="orange")
            if 150<=number[q][i][j]:
                lat_10m=int(min(lat))+ten_minute*j
                lon_10m=int(min(lon))+ten_minute*i
                ax.add_patch(patches.Rectangle((lon_10m, lat_10m),ten_minute,ten_minute,color='red'))
                p4=Rectangle((0, 0), 1, 1, fc="red")
    ax.legend([p1, p2,p3,p4], ['1~49','50~99','100~149','150~'],loc='upper left')
    ax.set_title(str(2006+q),fontsize=20)
plt.show()