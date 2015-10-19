#!/usr/bin/env /anaconda/bin/python
"""
Created on Fri May 29 09:46:15 2015
return forecast plot after choosing site and method on website
@author: zhaobin
"""
import cgi, cgitb
cgitb.enable()
import numpy as np 
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import netCDF4 
import datetime
import matplotlib.pyplot as plt
from MODULES import dm2dd,DIST,FIG_PLOT,whichArea,draw_basemap,stick_plot
from pylab import *
import matplotlib.tri as Tri
##################################################################
one_minute=1.0/60
sites=[[],[],[207,271],[247,278],[443,327],[264,369],[245,370]]   #different sites to plot
data=pd.read_csv('binned_td_FVCOM.csv')
temp_sf=pd.Series(data['MEAN_TEMP'])
time_sf=[]
for i in range(len(data)):
    time_sf.append(pd.Series(data['datet'])[i])
    pd.Series(data['datet'])[i]=datetime.datetime.strptime(pd.Series(data['datet'])[i],'%Y-%m-%d %H:%M:%S')
    time_sf[i]=datetime.datetime.strptime(time_sf[i],'%Y-%m-%d %H:%M:%S')  
time_sf=pd.Series(time_sf)
lat,lon=dm2dd(data['LATITUDE'],data['LONGITUDE']) #convert ddmm.m to dd.ddd
lat_i=np.arange(min(lat),max(lat),one_minute)
lon_i=np.arange(min(lon),max(lon),one_minute)
form = cgi.FieldStorage()
site = int(form.getvalue('site'))
method = form.getvalue('way')
lonsize=[lon_i[sites[site-1][0]]-5*one_minute,lon_i[sites[site-1][0]]+5*one_minute]
latsize=[lat_i[sites[site-1][1]]-5*one_minute,lat_i[sites[site-1][1]]+5*one_minute]
WAYS=[' ','Wind speed(m/s)','Wave height(m)','Bottom temperature(degC)','current']
way=['','wind','wave','bottom_temp','current']
if method=='wind':
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Forecasts/NECOFS_MET_FORECAST.nc'
    try:
        nc=netCDF4.Dataset(url)
        lons=nc.variables['XLONG'][:]
        lats=nc.variables['XLAT'][:] 
        time=nc.variables['Times'][:]
        u=nc.variables['U10'][:]
        v=nc.variables['V10'][:]
        Dist=[]
        for i in range(len(lons[:])):
            for j in range(len(lons[0])):
                Dist.append(DIST(lons[i][j],lats[i][j],lon_i[sites[site][0]],lat_i[sites[site][1]]))
        index_nearest=np.argmin(Dist)        # get smallest one
        index_one,index_two=index_nearest/len(lons[0]),index_nearest%(len(lons[0]))
        wind=[]
        for i in range(len(time)):
            U=u[i][index_one][index_two]
            V=v[i][index_one][index_two]
            s=np.sqrt(U*U+V*V)
            wind.append(s)
        T=[]
        for i in range(len(time)):
            t=''
            for j in range(len(time[i])):
                t=t+str(time[i][j])
            t=datetime.datetime.strptime(t,"%Y-%m-%d_%H:%M:%S")
            T.append(t)
        fig=FIG_PLOT(wind,T,site,method,WAYS)
        plt.savefig('/var/www/html/ioos/sf/fig/site'+str(site)+'_'+method+'.png')
    except:
        nc1='Model doesn`t work now.'
if method=='wave':#1 means wind,2 means wave,3 means bottom temp,4 means current,5 means all
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Forecasts/NECOFS_WAVE_FORECAST.nc'
    try:
        nc=netCDF4.Dataset(url)
        lons=nc.variables['lon'][:]
        lats=nc.variables['lat'][:]
        time=nc.variables['time'][:]
        wave=nc.variables['hs'][:]  
        Dist=[]
        for i in range(len(lons[:])):
            Dist.append(DIST(lons[i],lats[i],lon_i[sites[site-1][0]],lat_i[sites[site-1][1]]))
        index_nearest=np.argmin(Dist)
        T,WAVE=[],[]
        for i in range(len(time)):
            WAVE.append(wave[i][index_nearest])
            T.append(datetime.datetime(1858,11,17)+datetime.timedelta(days=int(time[i]))
                             +datetime.timedelta(seconds=int(time[i]%1*24*3600)))   
        fig=FIG_PLOT(WAVE,T,site,method,WAYS)
        plt.savefig('/var/www/html/ioos/sf/fig/site'+str(site)+'_'+method+'.png')
    except:
        nc1='Model doesn`t work now.'
if method=='bottom temperature':
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'
    try:
        nc=netCDF4.Dataset(url)
        lons=nc.variables['lon'][:]
        lats=nc.variables['lat'][:]
        time=nc.variables['time'][:]
        temp=nc.variables['temp']   
        Dist=[]
        for i in range(len(lons[:])):
            Dist.append(DIST(lons[i],lats[i],lon_i[sites[site-1][0]],lat_i[sites[site-1][1]]))
        index_nearest=np.argmin(Dist)
        T,TEMP_f=[],[]
        TEMP=temp[:,39,index_nearest]
        for i in range(len(time)):   
            T.append(datetime.datetime(1858,11,17)+datetime.timedelta(days=int(time[i]))
                +datetime.timedelta(seconds=int(time[i]%1*24*3600)))
        for i in range(len(TEMP)):
            TEMP_f.append(TEMP[i]*1.8+32)
        indx=[]
        for i in range(len(lon_i)):
            indx.append([])
            for j in range(len(lat_i)):
                indx[i].append([])
        for i in range(len(lat)): 
            n=whichArea(lat[i],lat_i)
            m=whichArea(lon[i],lon_i)
            indx[m][n].append(i)           # calculate index of data in 1 minute square
        now=datetime.datetime.now()+datetime.timedelta(hours=4) #utc to local time
        INDX=[]
        for i in range(10): 
            INDX.append([])   #2006~2015 is 10 years
        for i in indx[sites[site-1][0]][sites[site-1][1]]:
            time_sf[i]=time_sf[i].replace(year=now.year)  
            if T[0].date()<=time_sf[i].date()<now.date(): 
                for j in range(len(INDX)):
                    if pd.Series(data['datet'])[i].year==2006+j:
                        INDX[j].append(i)
        fig,ax,ax2=FIG_PLOT(TEMP,T,site,method,WAYS,TEMP_f)
        for j in range(len(INDX)):
            if len(INDX[j])>0:
                for i in INDX[j]:
                    temp_sf[i]=temp_sf[i]*1.8+32
                DATA={'time':time_sf[INDX[j]].values,'temp':temp_sf[INDX[j]].values}
                DATA=pd.DataFrame(DATA,index=time_sf[INDX[j]].values)
                DATA=DATA.sort_index()          
                ax2.plot(DATA['time'].values,DATA['temp'].values,label=str(2006+j),marker='o')
                ymin,ymax=ax.get_ylim()
                if (min(temp_sf[INDX[j]].values)-32)/1.8<ymin:
                    ax.set_ylim((min(temp_sf[INDX[j]].values)-32)/1.8,ymax)
                    ax2.set_ylim(min(temp_sf[INDX[j]].values),ymax*1.8+32) 
        plt.legend(loc='best')
        plt.savefig('/var/www/html/ioos/sf/fig/site'+str(site)+'_'+method+'.png')
    except:
        nc1='Model doesn`t work now.'
if method=='current':
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'
    try:
        nc=netCDF4.Dataset(url)
        lats = nc.variables['lat'][:]
        lons = nc.variables['lon'][:]
        latc = nc.variables['latc'][:]
        lonc = nc.variables['lonc'][:]
        time=nc.variables['time'][:]
        Dist=[]
        for i in range(len(lons[:])):
            Dist.append(DIST(lons[i],lats[i],lon_i[sites[site-1][0]],lat_i[sites[site-1][1]]))
        index_nearest=np.argmin(Dist)
        layer=0
        u = nc.variables['u'][:,layer,index_nearest]
        v = nc.variables['v'][:,layer,index_nearest]
        T=[]
        for i in range(len(time)):
            T.append(datetime.datetime(1858,11,17)+datetime.timedelta(days=int(time[i]))+datetime.timedelta(seconds=int(time[i]%1*24*3600)))
        q = stick_plot(T, u, v)
        maxvel =np.sqrt(np.max(abs(u))*np.max(abs(u))+np.max(abs(v))*np.max(abs(v)))
        qk = plt.quiverkey(q, 0.1, 0.85, maxvel,"%s N m$^{-2}$" % maxvel,labelpos='W', coordinates='axes')
        plt.xlabel('Time(month/day)',fontsize=20)    
        plt.title('Site'+str(site)+'current',fontsize=20)
        plt.savefig('/var/www/html/ioos/sf/fig/site'+str(site)+'_'+method+'.png')
    except:
        nc1='Model doesn`t work now.'
if method=='Current':
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'
    try:
        nc=netCDF4.Dataset(url)
        lons=nc.variables['lon'][:]
        lats=nc.variables['lat'][:]
        latc = nc.variables['latc'][:]
        lonc = nc.variables['lonc'][:]
        time=nc.variables['time'][:]  
        nv = nc.variables['nv'][:].T-1  # Get Connectivity array,T-1 means Transposed matrix
        h = nc.variables['h'][:] 
        start=datetime.datetime.now()+datetime.timedelta(hours=4)  #utc to local time
        T=[]
        for i in range(len(time)):
            T.append(datetime.datetime(1858,11,17)+datetime.timedelta(days=int(time[i]))
                   +datetime.timedelta(seconds=int(time[i]%1*24*3600)))
        itime=np.argmin(abs(np.array(T)-start))   #find nearest time
        daystr=T[itime]
        tri = Tri.Triangulation(lons,lats)
        ilayer = 0      #0 is surface,-1 is bottom
        u = nc.variables['u'][itime, ilayer, :]   
        v = nc.variables['v'][itime, ilayer, :]   
        ax= [lon_i[sites[site-1][0]]-5*one_minute,lon_i[sites[site-1][0]]+5*one_minute,
            lat_i[sites[site-1][1]]-5*one_minute,lat_i[sites[site-1][1]]+5*one_minute] # box you want to plot
        subsample=2
        ind = argwhere((lonc >= ax[0]) & (lonc <= ax[1]) & (latc >= ax[2]) & (latc <= ax[3]))  # find velocity points in bounding box
        np.random.shuffle(ind)
        Nvec = int(len(ind) / subsample)
        idv = ind[:Nvec]
        indx=argwhere((lons >= ax[0]) & (lons <= ax[1]) & (lats >= ax[2]) & (lats <= ax[3]))  # find  points in bounding box
        #plot picture
        fig=plt.figure(figsize=(12,10))
        ax=fig.add_subplot(111,aspect=(1.0/cos(mean(lat)*pi/180.0)))
        draw_basemap(fig, ax, lonsize, latsize,interval_lon=5*one_minute, interval_lat=5*one_minute)
        if int(-min(h[indx]))>-5:
            levels=arange(int(-max(h[indx]))-5,0,1)
        else:
            levels=arange(int(-max(h[indx]))-5,int(-min(h[indx]))+1,1)
        tricontourf(tri,-h,levels=levels,shading='faceted',cmap=plt.cm.gist_earth)
        gca().patch.set_facecolor('0.5')
        cbar=colorbar()
        cbar.set_label('Water Depth (m)')
        Q = quiver(lonc[idv],latc[idv],u[idv],v[idv],scale=10)
        maxvel =np.sqrt(np.max(abs(u[idv]))*np.max(abs(u[idv]))+np.max(abs(v[idv]))*np.max(abs(v[idv])))
        maxstr='%3.1f m/s' % maxvel
        qk = quiverkey(Q,0.92,0.08,maxvel,maxstr,labelpos='W')
        title('Site '+str(site)+' , surface ,'+str(daystr-datetime.timedelta(hours=4))[0:-3]+'(local time)',fontsize=15)
        plt.savefig('/var/www/html/ioos/sf/fig/site'+str(site)+'_'+method+'.png')
    except:
        nc1='Model doesn`t work now.'
if method=='Bottom temp':
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'
    try:
        nc=netCDF4.Dataset(url)
        lons=nc.variables['lon'][:]
        lats=nc.variables['lat'][:]
        time=nc.variables['time'][:]
        temp=nc.variables['temp']   
        start=datetime.datetime.now()+datetime.timedelta(hours=4)
        T=[]
        for i in range(len(time)):
            T.append(datetime.datetime(1858,11,17)+datetime.timedelta(days=int(time[i]))
                 +datetime.timedelta(seconds=int(time[i]%1*24*3600)))
        itime=np.argmin(abs(np.array(T)-start))   #find nearest time
        daystr=T[itime]
        BIN= [lon_i[sites[site-1][0]]-7*one_minute,lon_i[sites[site-1][0]]+7*one_minute,
                lat_i[sites[site-1][1]]-7*one_minute,lat_i[sites[site-1][1]]+7*one_minute] # box you want to plot,here 7>5 want to remove boundaries in ax.
        indx=argwhere((lons >= BIN[0]) & (lons <= BIN[1]) & (lats >= BIN[2]) & (lats <= BIN[3]))  # find  points in bounding box
        INDX=[]    
        for i in range(len(indx)):
            INDX.append(indx[i][0])
        TEMP=temp[itime,-1,INDX]   
        tri = Tri.Triangulation(lons[INDX],lats[INDX])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        draw_basemap(fig, ax, lonsize, latsize,interval_lon=5*one_minute, interval_lat=5*one_minute)
        levels=np.arange(round(min(TEMP),1),round(max(TEMP),1),0.1)
        tricontourf(tri,TEMP,levels=levels,cmap=plt.cm.rainbow)
        gca().patch.set_facecolor('0.5')
        cbar=colorbar()
        cbar.set_label('Temperature(degC)')
        plt.title('Site'+str(site)+',bottom temperature,'+str(daystr-datetime.timedelta(hours=4))[0:-3]+'(local time)',fontsize=20)
        plt.savefig('/var/www/html/ioos/sf/fig/site'+str(site)+'_'+method+'.png')
    except:
        nc1='Model doesn`t work now.'
if method=='Wave-height':
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Forecasts/NECOFS_WAVE_FORECAST.nc'
    try:
        nc=netCDF4.Dataset(url)
        lons=nc.variables['lon'][:]
        lats=nc.variables['lat'][:]
        time=nc.variables['time'][:]
        wave=nc.variables['hs'][:]
        start=datetime.datetime.now()+datetime.timedelta(hours=4)
        T=[]
        for i in range(len(time)):
            T.append(datetime.datetime(1858,11,17)+datetime.timedelta(days=int(time[i]))
                 +datetime.timedelta(seconds=int(time[i]%1*24*3600)))
        itime=np.argmin(abs(np.array(T)-start))   #find nearest time
        daystr=T[itime]
        BIN= [lon_i[sites[site-1][0]]-7*one_minute,lon_i[sites[site-1][0]]+7*one_minute,
                lat_i[sites[site-1][1]]-7*one_minute,lat_i[sites[site-1][1]]+7*one_minute] # box you want to plot,here 7>5 want to remove boundaries in ax.
        indx=argwhere((lons >= BIN[0]) & (lons <= BIN[1]) & (lats >= BIN[2]) & (lats <= BIN[3]))  # find  points in bounding box
        INDX=[]    
        for i in range(len(indx)):
            INDX.append(indx[i][0])
        WAVE=wave[itime,INDX]   
        tri = Tri.Triangulation(lons[INDX],lats[INDX])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        draw_basemap(fig, ax, lonsize, latsize,interval_lon=5*one_minute, interval_lat=5*one_minute)
        levels=np.arange(round(min(WAVE),3),round(max(WAVE),3),0.001)
        tricontourf(tri,WAVE,levels=levels,cmap=plt.cm.rainbow)
        gca().patch.set_facecolor('0.5')
        cbar=colorbar()
        cbar.set_label('Wave height(m)')
        plt.title('Site'+str(site)+',wave height,'+str(daystr-datetime.timedelta(hours=4))[0:-3]+'(local time)',fontsize=20)
        #plt.savefig('/var/www/html/ioos/sf/fig/site'+str(site)+'_'+method+'.png')
    except:
        nc1='Model doesn`t work now.'
try:
    if nc is not False:
        print "Content-type:text/html\r\n\r\n"
        print "<!DOCTYPE html PUBLIC '-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd'>"
        print "<html>"
        print "<head>"
        print "<title>Forecast</title>"
        print "</head>"
        print "<body>" 
        print "<img src='http://comet.nefsc.noaa.gov/ioos/sf/fig/site%s_%s.png' width='800' height='800' /> " %(site,method)
        print "</body>"
        print "</html>"
except:
    print "Content-type:text/html\r\n\r\n"
    print "<!DOCTYPE html PUBLIC '-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd'>"
    print "<html>"
    print "<head>"
    print "<title>Forecast</title>"
    print "</head>"
    print "<body>" 
    print "<font size='6' color='#FF0000'>Model can`t use now</font>"
    print "</body>"
    print "</html>"

