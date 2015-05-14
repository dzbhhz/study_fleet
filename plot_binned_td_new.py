# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:56:21 2013

@author: jmanning and zhaobin
"""

from pylab import *
import matplotlib.pyplot as plt
from pandas import *
import numpy as np
from datetime import datetime
from conversions_old import dm2dd #/home3/ocn/jmanning/py/conversions.py
from mpl_toolkits.basemap import Basemap
import os
import sys


# read in data by first using "/home3/ocn/jmanning/lob/r_and_d/hoey/write_binned.sql
fn='binned_td_2015_06.csv' # where I have added a header line after exporting with sqldump -u jmanning/emolt#00 -o /data5/jmanning/hoey/binned_td_2013_08 -f  /home3/ocn/jmanning/lob/r_and_d/hoey/write_binned.sql -D','
direct='/net/data5/jmanning/hoey/' # this is where I stored this imput file
PATH_IMG = 'animations/' # location of animations
def parse(datet):
   #print datet[0:10],datet[11:13],datet[14:16]
    dt=datetime.strptime(datet,'%Y-%m-%d:%H:%M')
    return dt
df=read_csv(direct+fn,sep=',',skiprows=0,parse_dates={'datet':[2]},index_col='datet',date_parser=parse,names=['LATITUDE','LONGITUDE','ROUND_DATE_TIME','OBSERVATIONS_TEMP','MEAN_TEMP','MIN_TEMP','MAX_TEMP','STD_DEV_TEMP','OBSERVATIONS_DEPTH','MEAN_DEPTH','MIN_DEPTH','MAX_DEPTH','STD_DEV_DEPTH','nan'])
for i in range(len(df)):
    if df['MEAN_DEPTH'][i]=='                                                                                                                               ':
        df['MEAN_DEPTH'][i]=0
    else:
        df['MEAN_DEPTH'][i]=float(df['MEAN_DEPTH'][i])
    if df['MIN_DEPTH'][i]=='                                                                                                                               ':
        df['MIN_DEPTH'][i]=df['MEAN_DEPTH'][i]
    else:
        df['MIN_DEPTH'][i]=float(df['MIN_DEPTH'][i])
    if df['MAX_DEPTH'][i]=='                                                                                                                               ':
        df['MAX_DEPTH'][i]=df['MEAN_DEPTH'][i]
    else:
        df['MAX_DEPTH'][i]=float(df['MAX_DEPTH'][i])
#convert ddmm.m to dd.ddd
[la,lo]=dm2dd(df['LATITUDE'],df['LONGITUDE'])

# Now we have a data frame, we can start plotting it
# make basemap
llLat=min(la)
urLat=max(la)
llLon=min(lo)
urLon=max(lo)

# Summary plot

m = Basemap(projection='cyl', llcrnrlat=llLat, urcrnrlat=urLat, \
                    llcrnrlon=llLon, urcrnrlon=urLon,resolution='h')
m.drawcoastlines()
m.fillcontinents(color='gray')
m.drawparallels(np.arange(round(llLat), round(urLat),1.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(round(llLon), round(urLon),1.), labels=[0, 0, 0, 1])

plot(lo,la,'r.')
title('bins with data 2006-2015')
plt.show()
savefig('/net/nwebserver/epd/ocean/MainPage/td/binned_td_hoey.png')

'''
# animation by month

df['month']=df.index.month
# Now for the animation by month
mth=df['month'][0] # intial month

fig=plt.figure()
m = Basemap(projection='cyl', llcrnrlat=llLat, urcrnrlat=urLat, \
                    llcrnrlon=llLon, urcrnrlon=urLon,resolution='h')
m.drawcoastlines()
m.fillcontinents(color='gray')
m.drawparallels(np.arange(round(llLat), round(urLat),1.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(round(llLon), round(urLon),1.), labels=[0, 0, 0, 1])
title(str(df.index.month[0])+'/'+str(df.index.year[0]))
for k in range(len(df)):
    if df['month'][k]==mth:
        plot(lo[k],la[k],'r.')
    else:
        mth=mth+1
        if mth==13:
            mth=1;
        plt.show()    
        savefig('/net/home3/ocn/jmanning/py/animations/'+str(df.index.year[k-1]*100+df.index.month[k-1])+'.png')    
        plt.clf()
        close('all')
        fig=plt.figure()
        m = Basemap(projection='cyl', llcrnrlat=llLat, urcrnrlat=urLat, \
                    llcrnrlon=llLon, urcrnrlon=urLon,resolution='h')
        m.drawcoastlines()
        m.fillcontinents(color='gray')
        m.drawparallels(np.arange(round(llLat), round(urLat),1.), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(round(llLon), round(urLon),1.), labels=[0, 0, 0, 1])
        title(str(df.index.month[k])+'/'+str(df.index.year[k]))
        print str(df.index.month[k])+'/'+str(df.index.year[k])
       
#make animation 
try:
    # Make .flv
    #anim_name = INPUT_FILENAME[0:-4] + '.flv' #Change extension to create other type of files; Make sure that ffmpeg supports it
    anim_name = '/net/nwebserver/epd/ocean/MainPage/td/'+fn[0:-4] + '.gif'
    #cmd = PATH_FFMPEG + 'ffmpeg.exe -i ' + PATH_IMG + '  -r 15 -b 614000 ' + PATH_ANIM + anim_name
    cmd = 'convert -delay 40 -loop 0 '+ PATH_IMG + '2*.png ' + anim_name
    os.system(cmd)
    print "Animation was created successfully."
 
except:
    print  "Could not create animation."
'''
'''
# Depth vs Time
plt.close('all')        
plot_date(df.index,-df['MEAN_DEPTH'])
ylim(-500,0)
title('Depth of observations vs time')
savefig('/net/nwebserver/epd/ocean/MainPage/td/binned_td_hoey_depth.png')

#expoert data to file
df.to_csv('/net/nwebserver/epd/ocean/MainPage/td/binned_td_hoey.csv',float_format='%10.2f')

#distribution plots
f,ax1 = plt.subplots(2)#sharex=True
bins=range(3,22)
ax1[0].hist(df['MEAN_TEMP'],bins)
ax1[0].set_title('Distributions')
ax1[0].set_ylabel('# of Obs')
ax1[0].set_xlabel('temperature (degC)')
bins=range(0,320,20)
ax1[1].hist(df['MEAN_DEPTH'],bins)
ax1[1].set_xlabel('Depths of Observations')
ax1[1].set_ylabel('# of Obs')
'''
plt.show()
