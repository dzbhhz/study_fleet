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
import matplotlib.pyplot as plt
from MODULES import DIST,FIG_PLOT_emolt_30yr,get_emolt_data
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
dist_emolt,data_emolt,emolt_name=get_emolt_data(alon,alat)
start_time,end_time=min(data_emolt['time']),max(data_emolt['time'])
if method=='bottom temperature':
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3/mean'
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
        t1=(start_time-datetime.datetime(1858,11,17)).total_seconds()/3600/24
        t2=(end_time-datetime.datetime(1858,11,17)).total_seconds()/3600/24
        nearest_s,nearest_e=np.argmin(abs(time-t1)),np.argmin(abs(time-t2))
        T=[]    
        for i in time[nearest_s:nearest_e]:
            t=datetime.datetime(1858,11,17)+datetime.timedelta(days=int(i))
            T.append(t.date())   
        TEMP=temp[nearest_s:nearest_e,-1,index_nearest]
        TEMP_f=[]
        for i in range(len(TEMP)):
            TEMP_f.append(TEMP[i]*1.8+32)
        fig,ax,ax2=FIG_PLOT_emolt_30yr(TEMP,T,alon,alat,method,WAYS,TEMP_f)
        ax.plot(data_emolt['time'],data_emolt['temperature'])
        ax2.plot(data_emolt['time'],data_emolt['temperature']*1.8+32,label='eMOLT')
        xmin,xmax=ax.get_ylim()
        ax2.set_ylim(xmin*1.8+32,xmax*1.8+32)
        plt.legend(loc='best') 
        plt.savefig('/var/www/html/ioos/sf/fig/'+method+'.png')
    except:
        nc1='Model does`t work now'
try:
    if nc is not False:
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
