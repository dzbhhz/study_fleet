#!/usr/bin/env /anaconda/bin/python
"""
Created on Fri May 29 09:46:15 2015

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
import pytz
import matplotlib.pyplot as plt
from MODULES import dm2dd,DIST,FIG_PLOT_emolt,whichArea,draw_basemap,stick_plot,get_emolt_data,getobs_tempsalt,gmt_to_loc
from pylab import *
import matplotlib.tri as Tri
from matplotlib.mlab import griddata
##################################################################
one_minute=1.0/60
form = cgi.FieldStorage()
alat = float(form.getvalue('alat'))
alon = float(form.getvalue('alon'))
method = form.getvalue('way')
lonsize=[alon-5*one_minute,alon+5*one_minute]
latsize=[alat-5*one_minute,alat+5*one_minute]
WAYS=[' ','Wind speed(m/s)','Wave height(m)','Bottom temperature(degC)','current']
way=['','wind','wave','bottom_temp','current']
if method=='bottom temperature':
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc'
    try:
        nc=netCDF4.Dataset(url)
        lons=nc.variables['lon'][:]
        lats=nc.variables['lat'][:]
        time=nc.variables['time'][:]
        temp=nc.variables['temp']   
        Dist=[]
        for i in range(len(lons[:])):
            Dist.append(DIST(lons[i],lats[i],alon,alat))
        index_nearest=np.argmin(Dist)
        depth_fvcom=nc.variables['h'][index_nearest]  #get depth of fvcom
        T,TEMP_f=[],[]
        TEMP=temp[:,-1,index_nearest]
        for i in range(len(time)):   
            t=datetime.datetime(1858,11,17)+datetime.timedelta(days=int(time[i]))+datetime.timedelta(seconds=int(time[i]%1*24*3600))
            t=gmt_to_loc(t)
            T.append(t)   
        for i in range(len(T)):
            T[i]=T[i].replace(year=2000,tzinfo=None)
        for i in range(len(TEMP)):
            TEMP_f.append(TEMP[i]*1.8+32)
        fig,ax,ax2=FIG_PLOT_emolt(TEMP,T,alon,alat,method,WAYS,TEMP_f)
        dist_emolt,data_emolt,emolt_name=get_emolt_data(alon,alat)
        i_emolt=[] 
        start,end=min(data_emolt['time']).year,max(data_emolt['time']).year
        dif_year=end-start+1
        for i in range(dif_year):
            i_emolt.append([]) 
            c=data_emolt[data_emolt.time>=datetime.datetime(start+i,1,1,0,0,0)]
            d=c[c.time<=datetime.datetime(start+i,12,31,23,59,59)].index.values
            i_emolt[i].extend(d)
        I_emolt,INDX,TIME,TEMP=[],[],[],[]
        for i in range(dif_year):
            I_emolt.append([])
            for j in range(len(i_emolt[i])):
                TIME.append(data_emolt['time'][i_emolt[i][j]].replace(year=2000))
                TEMP.append(data_emolt['temperature'][i_emolt[i][j]]*1.8+32)
                I_emolt[i].append(i_emolt[i][j])
                INDX.append(i_emolt[i][j])
        for i in range(dif_year):
            for j in range(len(i_emolt[i])):
                if T[0]>TIME[i_emolt[i][j]] or TIME[i_emolt[i][j]]>T[-1]:
                    I_emolt[i].remove(i_emolt[i][j])
                    INDX.remove(i_emolt[i][j])
        TEMP,TIME=pd.Series(TEMP),pd.Series(TIME)
        emolt_n=len(TEMP[INDX])
        if emolt_n>0:
            for i in range(dif_year):
                if len(TEMP[I_emolt[i]])>1:
                    data={'temp':TEMP[I_emolt[i]],'time':TIME[I_emolt[i]]}
                    Data=pd.DataFrame(data)
                    ax2.plot(Data['time'].values,Data['temp'].values,label=str(start+i))
            ymin,ymax=ax.get_ylim()
            if min(TEMP[INDX].values)<ymin*1.8+32:
                ax.set_ylim((min(TEMP[INDX].values)-32)/1.8,ymax)
                ax2.set_ylim(min(TEMP[INDX].values),ymax*1.8+32)
            ymin,ymax=ax.get_ylim()
            if max(TEMP[INDX].values)>ymax*1.8+32:
                ax.set_ylim(ymin,(max(TEMP[INDX].values)-32)/1.8)
                ax2.set_ylim(ymin*1.8+32,max(TEMP[INDX].values))
        xmin,xmax=ax.get_xlim()
        ax.set_xlim(xmin-1,xmax)
        ax2.set_xlim(xmin-1,xmax)
        plt.legend(loc='best') 
        plt.savefig('/var/www/html/ioos/sf/fig/'+method+'.png')
    except:
        nc1='Model does`t work now'
if method=='current':
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc'
    try:
        nc=netCDF4.Dataset(url)
        lats = nc.variables['lat'][:]
        lons = nc.variables['lon'][:]
        latc = nc.variables['latc'][:]
        lonc = nc.variables['lonc'][:]
        time=nc.variables['time'][:]
        Dist=[]
        for i in range(len(lons[:])):
            Dist.append(DIST(lons[i],lats[i],alon,alat))
        index_nearest=np.argmin(Dist)
        layer=0
        u = nc.variables['u'][:,layer,index_nearest]
        v = nc.variables['v'][:,layer,index_nearest]
        T=[]
        for i in range(len(time)):
            t=datetime.datetime(1858,11,17)+datetime.timedelta(days=int(time[i]))+datetime.timedelta(seconds=int(time[i]%1*24*3600))
            t=gmt_to_loc(t)
            T.append(t)
        q = stick_plot(T, u, v)
        maxvel =np.sqrt(np.max(abs(u))*np.max(abs(u))+np.max(abs(v))*np.max(abs(v)))
        qk = plt.quiverkey(q, 0.1, 0.85, maxvel,"%s m/s$^{-2}$" % round(maxvel,2),labelpos='W', coordinates='axes')
        plt.xlabel('Time(month/day)',fontsize=20)    
        plt.title('Loc:'+str(round(alon,1))+','+str(round(alat,1))+',current',fontsize=20)
        plt.savefig('/var/www/html/ioos/sf/fig/'+method+'.png')
    except:
        nc1='Model does`t work now'
if method=='Current':
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc'
    try:
        nc=netCDF4.Dataset(url)
        lons=nc.variables['lon'][:]
        lats=nc.variables['lat'][:]
        latc = nc.variables['latc'][:]
        lonc = nc.variables['lonc'][:]
        time=nc.variables['time'][:]  
        h = nc.variables['h'][:] 
        start=datetime.datetime.now(pytz.timezone("America/New_York"))
        T=[]
        for i in range(len(time)):
            t=datetime.datetime(1858,11,17)+datetime.timedelta(days=int(time[i]))+datetime.timedelta(seconds=int(time[i]%1*24*3600))
            t=gmt_to_loc(t)
            T.append(t)
        itime=np.argmin(abs(np.array(T)-start))   #find nearest time
        daystr=T[itime]
        tri = Tri.Triangulation(lons,lats)
        ilayer = 0      #0 is surface,-1 is bottom
        u = nc.variables['u'][itime, ilayer, :]   
        v = nc.variables['v'][itime, ilayer, :]   
        CX= [alon-7*one_minute,alon+7*one_minute,alat-7*one_minute,alat+7*one_minute] # box you want to plot
        subsample=2
        ind = argwhere((lonc >= CX[0]) & (lonc <= CX[1]) & (latc >= CX[2]) & (latc <= CX[3]))  # find velocity points in bounding box
        np.random.shuffle(ind)
        Nvec = int(len(ind) / subsample)
        idv = ind[:Nvec]
        indx=argwhere((lons >= CX[0]) & (lons <= CX[1]) & (lats >= CX[2]) & (lats <= CX[3]))  # find  points in bounding box
        #plot picture
        fig=plt.figure(figsize=(12,10))
        ax=fig.add_subplot(111)  #,aspect=(1.0/cos(mean(lat)*pi/180.0)))
        draw_basemap(fig, ax, lonsize, latsize,interval_lon=5*one_minute, interval_lat=5*one_minute)
        diff=int((max(h[indx])-min(h[indx]))/3)
        levels=arange(int(-max(h[indx]))-diff,int(-min(h[indx]))+diff,1)
        tricontourf(tri,-h,levels=levels,shading='faceted',cmap=plt.cm.gist_earth)
        gca().patch.set_facecolor('0.5')
        cbar=colorbar()
        cbar.set_label('Water Depth (m)')
        Q = quiver(lonc[idv],latc[idv],u[idv],v[idv],scale=10)
        maxvel =np.sqrt(np.max(abs(u[idv]))*np.max(abs(u[idv]))+np.max(abs(v[idv]))*np.max(abs(v[idv])))
        maxstr='%3.1f m/s' % maxvel
        qk = quiverkey(Q,0.92,0.08,maxvel,maxstr,labelpos='W')
        title('surface ,'+str(daystr)[0:-6]+'(local time)',fontsize=15)
        plt.savefig('/var/www/html/ioos/sf/fig/'+method+'.png')
    except:
        nc1='Model does`t work now'
if method=='Bottom temp':
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'
    try:
        nc=netCDF4.Dataset(url)
        lons=nc.variables['lon'][:]
        lats=nc.variables['lat'][:]
        time=nc.variables['time'][:]
        depth=-nc.variables['h'][:]
        temp=nc.variables['temp']   
        start=datetime.datetime.now(pytz.timezone("America/New_York"))
        T=[]
        for i in range(len(time)):
            t=datetime.datetime(1858,11,17)+datetime.timedelta(days=int(time[i]))+datetime.timedelta(seconds=int(time[i]%1*24*3600))
            t=gmt_to_loc(t)
            T.append(t)
        itime=np.argmin(abs(np.array(T)-start))   #find nearest time
        daystr=T[itime]
        BIN= [alon-20*one_minute,alon+20*one_minute,alat-20*one_minute,alat+20*one_minute] # box you want to plot,here 20>5 want to remove boundaries in ax.
        indx=argwhere((lons >= BIN[0]) & (lons <= BIN[1]) & (lats >= BIN[2]) & (lats <= BIN[3]))  # find  points in bounding box
        INDX=[]    
        for i in range(len(indx)):
            INDX.append(indx[i][0])
        TEMP=temp[itime,-1,INDX]
        H=depth[INDX]   
        tri = Tri.Triangulation(lons[INDX],lats[INDX])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        draw_basemap(fig, ax, lonsize, latsize,interval_lon=5*one_minute, interval_lat=5*one_minute)
        Bin_inside= [alon-7*one_minute,alon+7*one_minute,alat-7*one_minute,alat+7*one_minute]
        indx_inside=argwhere((lons >= Bin_inside[0]) & (lons <= Bin_inside[1]) & (lats >= Bin_inside[2]) & (lats <= Bin_inside[3]))  # find  points in bounding box
        INDX_inside=[]    
        for i in range(len(indx_inside)):
            INDX_inside.append(indx_inside[i][0])
        TEMP_inside=temp[itime,-1,INDX_inside]
        H_inside=depth[INDX_inside]
        levels=np.arange(round(min(TEMP_inside),2),round(max(TEMP_inside),2),0.01)
        h_start,h_end=round(min(H_inside)),round(max(H_inside))
        levels_h=np.arange(h_start,h_end,(h_end-h_start)/5)
        cs=tricontourf(tri,TEMP,levels=levels,cmap=plt.cm.rainbow)
        cs_line=tricontour(tri,H,levels=levels_h,linewidths=0.5, colors='k')
        plt.clabel(cs_line, cs_line.levels, fontsize=10)
        gca().patch.set_facecolor('0.5')
        cbar=colorbar(cs)
        cbar.set_label('Temperature(degC)')
        plt.title('Bottom temperature,'+str(daystr)[0:-6]+'(local time)',fontsize=15)
        plt.savefig('/var/www/html/ioos/sf/fig/'+method+'.png')
    except:
        nc1='Model does`t work now'
try:
    if type(nc) is netCDF4._netCDF4.Dataset:
        print "Content-type:text/html\r\n\r\n"
        print "<!DOCTYPE html PUBLIC '-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd'>"
        print "<html>"
        print "<head>"
        print "<title>Forecast</title>"
        print "</head>"
        print "<body>" 
        if method=='bottom temperature':
            if dist_emolt>20:
                print '<font size="6" color="#FF0000">WARNING:Site selected is %s km from nearest eMOLT site(%s) </font>' %(round(dist_emolt,1),emolt_name)
                if emolt_n==0:
                    print "<p>"
                    print '<font size="6" color="#FF0000">WARNING:No emolt data in these six days</font>' 
            else:
                print 'Site selected is %s km from nearest eMOLT site(%s)' %(round(dist_emolt,1),emolt_name)
                if emolt_n==0:
                    print "<p>"
                    print '<font size="6" color="#FF0000">WARNING:No emolt data in these six days</font>'
            print "<p>"
            print 'Depth of '+str(emolt_name)+' is '+str(round(data_emolt['depth'][0],1))+'m'
            print "<p>"        
            print 'Depth of selected site is '+str(round(depth_fvcom,1))+'m'
            print "<p>"
        print "<img src='http://comet.nefsc.noaa.gov/ioos/sf/fig/%s.png' width='800' height='800' /> " %method
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
