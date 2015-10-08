# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 09:15:48 2015

@author: zhaobin
"""
'plot Modeled - observed temperature in places(observations>limit) '
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
import numpy as np
import pylab  
from conversions_old import dm2dd
from turtleModule import draw_basemap,whichArea
########################################################
one_minute=1.0/60
limit=100                   #set a limit you want to plot
sites=[[185,281],[193,273],[207,271],[247,278],[443,327],[264,369],[245,370]]   #different sites to plot
data=pd.read_csv('binned_td_FVCOM.csv')
for i in range(len(data)):
    data['datet'][i]=datetime.strptime(data['datet'][i], "%Y-%m-%d %H:%M:%S")
lat,lon=dm2dd(data['LATITUDE'],data['LONGITUDE']) #convert ddmm.m to dd.ddd
if limit==10:
    lonsize=[min(lon),max(lon)]
    latsize=[min(lat),max(lat)]
if limit==100:    #cut range of pic
    lonsize=[-72,-67]
    latsize=[40,43]
lat_i=np.arange(min(lat),max(lat),one_minute)   
lon_i=np.arange(min(lon),max(lon),one_minute)    #1.0/60 is 1 minute
number=[]
for i in range(len(lon_i)):
    j=[0]*len(lat_i)
    number.append(j)
index=[]
for i in range(len(lon_i)):
    index.append([])
    for j in range(len(lat_i)):
        index[i].append([])
for i in range(len(lat)):
    n=whichArea(lat[i],lat_i)
    m=whichArea(lon[i],lon_i)
    number[m][n]+=1                                      # calculate number in 1 minute bin
    index[m][n].append(i)                         #calculate index in 1-minute bin
Number=[]
for i in range(len(lat)):
    n=whichArea(lat[i],lat_i)
    m=whichArea(lon[i],lon_i)
    if number[m][n]>limit:                         #number>limit in 1 minute bin
        Number.append([m,n])
unique_Number = [list(x) for x in set(tuple(x) for x in Number)]   #get unique list
temp_mod=[]
temp_obs=[]
for i in range(len(lon_i)):
    temp_mod.append([])
    for j in range(len(lat_i)):
        temp_mod[i].append([])
for i in range(len(lon_i)):
    temp_obs.append([])
    for j in range(len(lat_i)):
        temp_obs[i].append([])
for i in range(len(temp_mod)):
    for j in range(len(temp_mod[i])):
        for q in index[i][j]:
            temp_mod[i][j].append(data['fvcom_temp'][q])
            temp_obs[i][j].append(data['MEAN_TEMP'][q])
Temp_mod,Temp_obs=[],[]
for i in range(len(lon_i)):
    Temp_mod.append([])
    for j in range(len(lat_i)):
        Temp_mod[i].append([])
for i in range(len(lon_i)):
    Temp_obs.append([])
    for j in range(len(lat_i)):
        Temp_obs[i].append([])
for i in range(len(Temp_mod)):
    for j in range(len(Temp_mod[i])):
        Temp_mod[i][j].append(np.mean(temp_mod[i][j]))
        Temp_obs[i][j].append(np.mean(temp_obs[i][j]))

fig=plt.figure(figsize=(12,10))
ax=fig.add_subplot(111)
draw_basemap(fig,ax,lonsize,latsize)
for i in range(len(unique_Number)):   
    diff=Temp_mod[unique_Number[i][0]][unique_Number[i][1]]-Temp_obs[unique_Number[i][0]][unique_Number[i][1]][0]
    if diff>0:
        plt.scatter(lon_i[unique_Number[i][0]],lat_i[unique_Number[i][1]],diff*50,color='r',edgecolors='black')
    if diff<0:
        plt.scatter(lon_i[unique_Number[i][0]],lat_i[unique_Number[i][1]],abs(diff*50),color='b',edgecolors='black')  #    *50 want to plot dots bigger
l1=plt.scatter([],[],s=25,edgecolors='none',color='r')
l2=plt.scatter([],[],s=50,edgecolors='none',color='r')
l3=plt.scatter([],[],s=100,edgecolors='none',color='r')
l4=plt.scatter([],[],s=150,edgecolors='none',color='r')
labels = ["0.5", "1.0", "2.0", "3.0"]
leg = plt.legend([l1, l2, l3, l4], labels,loc='upper left',labelspacing=1,columnspacing=1,prop={'size':15},scatterpoints=1)
pylab.rcParams['legend.numpoints'] = 1
if limit==100:
    for i in range(len(sites)):
        plt.annotate(str(i), xy=(lon_i[sites[i][0]], lat_i[sites[i][1]]),xytext=(lon_i[sites[i][0]],lat_i[sites[i][1]]-0.3),
             arrowprops=dict(arrowstyle='->'),fontsize=12)
plt.title('Modeled - observed temperature in places(observations>'+str(limit)+')',fontsize=20)
plt.show()
