# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 15:56:28 2015
plot dots animated by month
@author: zhaobin
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation
from datetime import datetime,timedelta
import numpy as np
import sys
sys.path.append('../modules')
from conversions_old import dm2dd
from turtleModule import draw_basemap,dist,whichArea
#######################################################################################
data=pd.read_csv('binned_td_hoey.csv')
lat,lon=dm2dd(data['LATITUDE'],data['LONGITUDE'])  #convert ddmm.m to dd.ddd
lon=pd.Series(lon)
lat=pd.Series(lat)
lonsize=[min(lon),max(lon)]
latsize=[min(lat),max(lat)]
for i in data.index:
    data['datet'][i]=datetime.strptime(data['datet'][i],'%Y-%m-%d %H:%M:%S')
obsyears=[]
for i in range(10):
    obsyears.append([])
    for j in range(12):
        obsyears[i].append([])
for i in range(len(obsyears)):
    for j in range(len(data)):
        if data['datet'][j].year==2006+i:
            for q in range(12):
                if data['datet'][j].month==q+1:
                    obsyears[i][q].append(j)
obsYear=[]
for i in range(10):                     #10 is the number of obsYear 2006~2015
    for j in range(12):                #12 is number of month
        obsYear.append(obsyears[i][j])

k=[0]*12+[1]*12+[2]*12+[3]*12+[4]*12+[5]*12+[6]*12+[7]*12+[8]*12+[9]*12           #use for year
K=[1,2,3,4,5,6,7,8,9,10,11,12]*10                     #use for month

fig=plt.figure()
ax=fig.add_subplot(111)
draw_basemap(fig, ax, lonsize, latsize)

def animate(i):
    del ax.collections[:] 
    plt.scatter(np.array(lon[obsYear[i]]),np.array(lat[obsYear[i]]), marker='o',edgecolors='none',c='r', s=10, zorder=15)
    plt.title('year='+str(2006+k[i])+',month='+str(K[i])+'') 
anim = animation.FuncAnimation(fig, animate, frames=119, interval=1000)    
anim.save('animated_by_month1.mov', fps=2)
plt.show()

