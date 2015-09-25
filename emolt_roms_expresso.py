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
from dateutil import tz
import matplotlib.pyplot as plt
from MODULES import dm2dd,DIST,FIG_PLOT_emolt,whichArea,draw_basemap,stick_plot,get_emolt_data,getobs_tempsalt,rot2d,shrink,bbox2ij,gmt_to_loc
from pylab import *
from mpl_toolkits.basemap import Basemap
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
if method=='bottom temperature':
    url='http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/his_Best/ESPRESSO_Real-Time_v2_History_Best_Available_best.ncd'
    try:
        nc=netCDF4.Dataset(url)
        lons=nc.variables['lon_rho'][:]
        lats=nc.variables['lat_rho'][:]
        time=nc.variables['time'][:]
        temp=nc.variables['temp']
        Dist=[]
        for i in range(len(lons[:])):
            Dist.append(DIST(lons[i],lats[i],alon,alat))  #calculate distance
        index_nearest=np.argmin(Dist)        #find index of nearest distance
        index_one,index_two=index_nearest/len(lons[0]),index_nearest%len(lons[0])
        starttime=datetime.datetime.now(pytz.timezone("America/New_York"))
        T_i=[]
        for i in range(len(time)):
            t=datetime.datetime(2013,5,18,0,0,0)+datetime.timedelta(hours=int(time[i]))
            t=gmt_to_loc(t)
            T_i.append(t)
        itime=np.argmin(abs(np.array(T_i)-starttime))
        depth_model=nc.variables['h'][index_one][index_two]  #get depth of roms
        TEMP_f=[]
        TEMP=temp[itime-72:itime+72,0,index_one,index_two]   #get nearest roms temperature about recent 6 days.
        if np.isnan(TEMP).any()==False:
            T=T_i[itime-72:itime+72]
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
            if len(TEMP[INDX])>0:
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
            plt.savefig('/var/www/html/ioos/sf/fig/roms'+method+'.png')
    except:
        nc1='Model does`t work now'
if method=='current':
    url='http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/his_Best/ESPRESSO_Real-Time_v2_History_Best_Available_best.ncd'
    try:
        nc=netCDF4.Dataset(url)
        lons=nc.variables['lon_rho'][:]
        lats=nc.variables['lat_rho'][:]
        time=nc.variables['time'][:]  
        Dist=[]
        for i in range(len(lons[:])):
            Dist.append(DIST(lons[i],lats[i],alon,alat))  #calculate distance
        index_nearest=np.argmin(Dist)        #find index of nearest distance
        index_one,index_two=index_nearest/len(lons[0]),index_nearest%len(lons[0])
        layer=-1  #in ROMS,last layer is surface
        starttime=datetime.datetime.now(pytz.timezone("America/New_York"))
        T_i=[]
        for i in range(len(time)):
            t=datetime.datetime(2013,5,18,0,0,0)+datetime.timedelta(hours=int(time[i]))
            t=gmt_to_loc(t)
            T_i.append(t)
        itime=np.argmin(abs(np.array(T_i)-starttime))
        u = nc.variables['u'][itime-72:itime+72,layer,index_one,index_two]
        v = nc.variables['v'][itime-72:itime+72,layer,index_one,index_two]
        TIME=nc.variables['time'][itime-72:itime+72]
        T=T_i[itime-72:itime+72]
        q = stick_plot(T, u, v)
        maxvel =np.sqrt(np.max(abs(u))*np.max(abs(u))+np.max(abs(v))*np.max(abs(v)))
        qk = plt.quiverkey(q, 0.1, 0.85, maxvel,"%s m/s$^{-2}$" % round(maxvel,2),labelpos='W', coordinates='axes')
        plt.xlabel('Time(month/day)',fontsize=20)    
        plt.title('Location:'+str(alon)+','+str(alat)+',current',fontsize=20)
        plt.savefig('/var/www/html/ioos/sf/fig/roms'+method+'.png')
    except:
        nc1='Model does`t work now'
if method=='Current':
    url = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/his_Best/ESPRESSO_Real-Time_v2_History_Best_Available_best.ncd'
    try:
        nc = netCDF4.Dataset(url)
        time = nc.variables['time'][:]  #recent 6 days '
        lon_rho = nc.variables['lon_rho'][:]
        lat_rho = nc.variables['lat_rho'][:]
        bbox = [alon-one_minute*30, alon+one_minute*30, alat-one_minute*30, alat+one_minute*30]
        i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
        starttime=datetime.datetime.now(pytz.timezone("America/New_York"))
        T=[]
        for i in range(len(time)):
            t=datetime.datetime(2013,5,18,0,0,0)+datetime.timedelta(hours=int(time[i]))
            t=gmt_to_loc(t)
            T.append(t)
        tidx =np.argmin(abs(np.array(T)-starttime))   #find nearest time
        zlev = -1  # last layer is surface layer in ROMS
        u = nc.variables['u'][tidx, zlev, j0:j1, i0:(i1-1)]
        v = nc.variables['v'][tidx, zlev, j0:(j1-1), i0:i1]
        lon=lon_rho[(j0+1):(j1-1), (i0+1):(i1-1)]
        lat=lat_rho[(j0+1):(j1-1), (i0+1):(i1-1)]
        mask = 1 - nc.variables['mask_rho'][(j0+1):(j1-1), (i0+1):(i1-1)]
        ang = nc.variables['angle'][(j0+1):(j1-1), (i0+1):(i1-1)]
        u = shrink(u, mask.shape)
        v = shrink(v, mask.shape)  # average u,v to central rho points
        U, V = rot2d(u, v, ang)  # rotate grid_oriented u,v to east/west u,v
        basemap = Basemap(projection='merc',llcrnrlat=bbox[2]-0.3,urcrnrlat=bbox[3],llcrnrlon=bbox[0]-0.3,urcrnrlon=bbox[1]+0.3, lat_ts=30,resolution='i')
        fig1 = plt.figure(figsize=(15,10))
        ax = fig1.add_subplot(111)
        basemap.drawcoastlines()
        basemap.fillcontinents()
        basemap.drawcountries()
        basemap.drawstates()
        basemap.drawparallels(np.arange(int(bbox[2]-0.5),int(bbox[3]+0.5),0.5),labels=[1,0,0,0], linewidth=0)
        basemap.drawmeridians(np.arange(int(bbox[0]-0.5),int(bbox[1]+0.5),0.5),labels=[0,0,0,1], linewidth=0)
        x_rho, y_rho = basemap(lon,lat)
        spd = np.sqrt(U*U + V*V)
        nsub=3
        scale=0.05
        Q=basemap.quiver(x_rho[::nsub,::nsub],y_rho[::nsub,::nsub],U[::nsub,::nsub],V[::nsub,::nsub],scale=1.0/scale, zorder=1e35, width=0.005)
        maxvel=np.nanmax(spd)
        maxstr='%3.1f m/s' % maxvel
        qk = quiverkey(Q,0.92,0.08,maxvel,maxstr,labelpos='W')
        plt.title('Surface current,'+str(T[tidx])[0:-6]+'(local time)',fontsize=15)
        plt.savefig('/var/www/html/ioos/sf/fig/roms'+method+'.png')
    except:
        nc1='Model does`t work now'
if method=='Bottom temp':
    url = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/his_Best/ESPRESSO_Real-Time_v2_History_Best_Available_best.ncd'
    try:
        nc = netCDF4.Dataset(url)
        time = nc.variables['time'][:]  #recent 6 days '
        lon_rho = nc.variables['lon_rho'][:]
        lat_rho = nc.variables['lat_rho'][:]
        bbox = [alon-one_minute*10, alon+one_minute*10, alat-one_minute*10, alat+one_minute*10]
        i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
        starttime=datetime.datetime.now(pytz.timezone("America/New_York"))
        T=[]
        for i in range(len(time)):
            t=datetime.datetime(2013,5,18,0,0,0)+datetime.timedelta(hours=int(time[i]))
            t=gmt_to_loc(t)
            T.append(t)
        tidx =np.argmin(abs(np.array(T)-starttime))   #find nearest time
        zlev = 0  # first layer is bottom layer in ROMS
        temp = nc.variables['temp'][tidx, zlev, j0:j1, i0:i1]
        lon = nc.variables['lon_rho'][j0:j1,i0:i1]
        lat = nc.variables['lat_rho'][j0:j1,i0:i1]
        depth=-nc.variables['h'][j0:j1,i0:i1]
        lon_inside,lat_inside,TEMP,H=[],[],[],[]
        for i in range(len(temp)):
            for j in range(len(temp[i])):
                TEMP.append(temp[i][j])
                lon_inside.append(lon[i][j])
                lat_inside.append(lat[i][j])
                H.append(depth[i][j])
        tri = Tri.Triangulation(lon_inside,lat_inside)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        draw_basemap(fig, ax, lonsize,latsize,interval_lon=0.1, interval_lat=0.1)
        levels=(np.arange(int(min(TEMP)),int(max(TEMP))+1.1,0.2))
        h_start,h_end=round(min(H)),round(max(H))
        levels_h=np.arange(h_start,h_end,(h_end-h_start)/10)
        levels_h=[int(i) for i in levels_h]
        cs=tricontourf(tri,TEMP,levels=levels,cmap=plt.cm.rainbow)
        cs_line=tricontour(tri,H,levels=levels_h,linewidths=0.5, colors='k')
        plt.clabel(cs_line, cs_line.levels,inline=1,inline_spacing=0,fmt='%1.0f', fontsize=10)
        gca().patch.set_facecolor('0.5')
        cbar=colorbar(cs)
        cbar.set_label('Temperature(degC)')
        plt.title('Bottom temperature,'+str(T[tidx])[0:-6]+'(local time)',fontsize=20) 
        plt.savefig('/var/www/html/ioos/sf/fig/roms'+method+'.png')
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
