# -*- coding: utf-8 -*-
"""
Created on Thu May 21 10:51:29 2015
Plots both a seasonal mean cycle and individual years of Study Fleet data
@author: zhaobin
"""
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
import matplotlib as mpl
import math
import numpy as np
from conversions_old import dm2dd
from turtleModule import whichArea
def remove_seacycle(df,DF):
    qqq=[]
    for j in range(len(df)):
        for i in range(len(DF)):
            if DF.index[i].week==df.index[j].week:
                if DF['mean'][i]>0:
                    q=df['MEAN_TEMP'][j]-DF['mean'][i]
                else:                            #some time don`t have temperature
                    q=float('nan')
                qqq.append(q)  
    df['Mean']=pd.Series(qqq,index=df.index)
    return df 
########################################################################
one_minute=1.0/60   #1.0/60 is 1 minute
data=pd.read_csv('binned_td_FVCOM_ROMS.csv')
lat,lon=dm2dd(data['LATITUDE'],data['LONGITUDE']) #convert ddmm.m to dd.ddd
lat_i=np.arange(min(lat),max(lat),one_minute)   
lon_i=np.arange(min(lon),max(lon),one_minute)
for i in data.index:
    data['datet'][i]=datetime.strptime(data['datet'][i],'%Y-%m-%d %H:%M:%S') 
sites=[[185,281],[193,273],[207,271],[247,278],[443,327],[264,369],[245,370]]   #different sites to plot
for s in range(len(sites)):
    indx=[]
    for i in range(len(lat)):
        m=whichArea(lon[i],lon_i)
        n=whichArea(lat[i],lat_i)
        if m==sites[s][0] and n==sites[s][1]:
            indx.append(i)                     # calculate number in each minute bin
    # start here to calculate the mean for each day of the year
    daily_ave=[]
    dt=[]
    for i in range(365):
        daily_temps=[]
        for k in indx:
             if i==data['datet'][k].timetuple().tm_yday:   #yearday
                 daily_temps.append(data['MEAN_TEMP'][k])
        daily_ave.append(np.mean(daily_temps))
        dt.append(datetime(2000,1,1,0,0,0)+timedelta(days=i))
    value_day=pd.Series(daily_ave,index=dt)
    value_week=value_day.resample('W',how=['mean'],loffset=timedelta(hours=-12))
    indx_year=[]
    for j in range(10): #2006~2015 is 10
        indx_year.append([])
        for i in indx:  # for all cases of the data inside the bin
            if data['datet'][i].year==2006+j: 
                indx_year[j].append(i)       #get index of data in different years
    for i in indx: # for all cases of the data inside the bin
        for j in indx: 
            if data['datet'][i]==data['datet'][j] and i!=j:
                indx.remove(i)       #remove index of some data which are in same time and same site
                data['MEAN_TEMP'][j]=np.mean([data['MEAN_TEMP'][i],data['MEAN_TEMP'][i]])  #calculate average of some data which are in same time and same site
                for q in range(len(indx_year)):
                    for k in indx_year[q]:
                        if k==i:
                            indx_year[q].remove(i)     #remove index of some data which are in same time and same site                       
    values={'MEAN_TEMP':data['MEAN_TEMP'][indx].values}
    VALUES=pd.DataFrame(values,index=data['datet'][indx])
    value_without_season=remove_seacycle(VALUES,value_week)   #remove seasonal cycle
    fig=plt.figure(figsize=(15,10))
    ax=fig.add_subplot(211)
    for i in range(10):   #10 is 2006~2015
        if len(indx_year[i])>1:    #>1 want to plot a line not a dot
            obs=ax.plot(data['datet'][indx_year[i]].values,data['MEAN_TEMP'][indx_year[i]].values,c='black',linewidth=2)
            fvcom=ax.plot(data['datet'][indx_year[i]].values,data['fvcom_temp'][indx_year[i]].values,c='blue',linewidth=2,linestyle='--')
            if s<4:    #site1~4 in roms range
                roms=ax.plot(data['datet'][indx_year[i]].values,data['temp_roms'][indx_year[i]].values,c='red',linewidth=2,linestyle='--')
    ax.set_xticks([min(data['datet'][indx]),max(data['datet'][indx])])
    ax.legend(['observed','fvcom','roms'],loc='best', ncol=3, fancybox=True, shadow=True)
    ax.set_ylabel('Temperature(degC)',fontsize=15)
    ax.set_title('Temperature vs time',fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.grid(True)
    ax1=fig.add_subplot(212,sharex=ax)
    Indx=[]
    for i in range(10): 
        Indx.append([])
        for j in range(len(value_without_season)):
            if value_without_season.index[j].year==2006+i:
                Indx[i].append(j)
        if len(indx_year[i])>0:
            ax1.plot(data['datet'][indx_year[i]].values,value_without_season['Mean'][value_without_season.index[Indx[i]]].values,c='black')
    ax1.set_xlabel('Year',fontsize=20)
    ax1.set_ylabel('Temperature(degC)',fontsize=15)
    ax1.set_title('With seasonal cycle removed',fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax1.grid(True)
    plt.setp(ax.get_xticklabels(),visible=False)
    plt.savefig('remove_seasonal_'+str(s+1)+'.png')
    fig1=plt.figure(figsize=(15,10))
    ax2=fig1.add_subplot(111)
    ax2.plot(value_week.index,value_week.values,c='black',linestyle='--',linewidth=5,label='weekly ave')
    for i in range(10):
        for j in indx_year[i]:
            data['datet'][j]=data['datet'][j].replace(year=2000)
        if len(indx_year[i])>0:
            ax2.plot(data['datet'][indx_year[i]].values,data['MEAN_TEMP'][indx_year[i]].values,label=str(2006+i))
    plt.legend(loc='best')
    dates = mpl.dates.drange(value_day.index[min(np.where(value_day.notnull())[0])],
                    value_day.index[max(np.where(value_day.notnull())[0])], timedelta(days=30))
    dateFmt = mpl.dates.DateFormatter('%m')
    ax2.set_xticks(dates)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax2.xaxis.set_major_formatter(dateFmt)
    ax2.set_xlabel('Month',fontsize=20)
    ax2.set_ylabel('Temperature(degC)',fontsize=20)
    ax2.set_title('Temperature vs seasonal cycle',fontsize=30)
    ax2.grid(True)
    plt.savefig('seasonal_'+str(s+1)+'.png')
plt.show()
