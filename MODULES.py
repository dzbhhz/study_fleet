# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 13:47:53 2015

@author: zhaobin
"""
import netCDF4
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.dates as mdates
import datetime
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib import path 
from dateutil.parser import parse
from dateutil import tz
import pytz
import datetime as dt
def FIG_PLOT(wind,T,site,method,way,TEMP_f=0):
    #plot picture
    fig=plt.figure(figsize=(10,10))
    ax=fig.add_subplot(111)
    ax.plot(T,wind)
    dates = mpl.dates.drange(min(T),max(T),datetime.timedelta(days=1))
    dateFmt = mpl.dates.DateFormatter('%m/%d')
    ax.set_xticks(dates)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.xaxis.set_major_formatter(dateFmt)
    ax.xaxis.set_major_locator(mdates.DayLocator())
    plt.gcf().autofmt_xdate()
    ax.set_xlabel('Time(month/day)',fontsize=20)
    if method=='bottom temperature':
        ymin,ymax=ax.get_ylim()
        ax2=ax.twinx()
        ax2.plot(T,TEMP_f,label='model forecast')
        ax2.set_ylabel('Bottom temperature(degF)',fontsize=20)
        ax2.set_ylim(ymin*1.8+32,ymax*1.8+32)
        ax.set_xlabel('Time',fontsize=20)
    if method=="wind":
        ax.set_ylabel(method+'(m/s)',fontsize=20)
    if method=="wave":
        ax.set_ylabel(method+'(m)',fontsize=20)
    if method=="bottom temperature":
        ax.set_ylabel(method+'(degC)',fontsize=20)
    ax.set_title('Site_'+str(site),fontsize=30)
    if method=='bottom temperature':
        return fig,ax,ax2
    else:
        return fig
def FIG_PLOT_emolt(wind,T,alon,alat,method,way,TEMP_f=0):
    #plot picture
    fig=plt.figure(figsize=(10,10))
    ax=fig.add_subplot(111)
    ax.plot(T,wind,linewidth=5)
    dates = mpl.dates.drange(min(T),max(T),datetime.timedelta(days=1))
    dateFmt = mpl.dates.DateFormatter('%m/%d')
    ax.set_xticks(dates)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.xaxis.set_major_formatter(dateFmt)
    ax.xaxis.set_major_locator(mdates.DayLocator())
    plt.gcf().autofmt_xdate()
    ax.set_xlabel('Time(month/day)',fontsize=20)
    if method=='bottom temperature':
        ymin,ymax=ax.get_ylim()
        ax2=ax.twinx()
        ax2.plot(T,TEMP_f,label='model',linewidth=5)
        ax2.set_ylabel('Bottom temperature(degF)',fontsize=20)
        ax2.set_ylim(ymin*1.8+32,ymax*1.8+32)
        dates = mpl.dates.drange(min(T),max(T),datetime.timedelta(days=1))
        dateFmt = mpl.dates.DateFormatter('%m/%d')
        ax2.set_xticks(dates)
        ax2.tick_params(axis='both', which='major', labelsize=15)
        ax2.xaxis.set_major_formatter(dateFmt)
        ax2.xaxis.set_major_locator(mdates.DayLocator())
        plt.gcf().autofmt_xdate()
        ax2.set_xlabel('Time(month/day)',fontsize=20)
    if method=="wind":
        ax.set_ylabel(method+'(m/s)',fontsize=20)
    if method=="wave":
        ax.set_ylabel(method+'(m)',fontsize=20)
    if method=="bottom temperature":
        ax.set_ylabel(method+'(degC)',fontsize=20)
    ax.set_title('Loaction:'+str(alon)+','+str(alat),fontsize=30)
    if method=='bottom temperature':
        return fig,ax,ax2
    else:
        return fig
def FIG_PLOT_emolt_30yr(wind,T,alon,alat,method,way,TEMP_f=0):
    #plot picture
    fig=plt.figure(figsize=(10,10))
    ax=fig.add_subplot(111)
    ax.plot(T,wind,linewidth=5,label='MODEL')
    ax.set_xlabel('Time',fontsize=20)
    if method=='bottom temperature':
        ymin,ymax=ax.get_ylim()
        ax2=ax.twinx()
        ax2.plot(T,TEMP_f,label='model',linewidth=5)
        ax2.set_ylabel('Bottom temperature(degF)',fontsize=20)
        ax2.set_ylim(ymin*1.8+32,ymax*1.8+32)
        ax2.set_xlabel('Time',fontsize=20)
    if method=="wind":
        ax.set_ylabel(method+'(m/s)',fontsize=20)
    if method=="wave":
        ax.set_ylabel(method+'(m)',fontsize=20)
    if method=="bottom temperature":
        ax.set_ylabel(method+'(degC)',fontsize=20)
    ax.set_title('Loaction:'+str(alon)+','+str(alat),fontsize=30)
    if method=='bottom temperature':
        return fig,ax,ax2
    else:
        return fig
def angle_conversion(a):
    a = np.array(a)
    return a/180*np.pi
def DIST(lon1, lat1, lon2, lat2):
    R = 6371.004
    lon1, lat1 = angle_conversion(lon1), angle_conversion(lat1)
    lon2, lat2 = angle_conversion(lon2), angle_conversion(lat2)
    l = R*np.arccos(np.cos(lat1)*np.cos(lat2)*np.cos(lon1-lon2)+\
                        np.sin(lat1)*np.sin(lat2))
    return l
#convert degrees, minutes to decimal
def dm2dd(la,lo):

  def dm2dd_single(lat,lon):
    (a,b)=divmod(float(lat),100)   
    aa=int(a)
    bb=float(b)
    lat_value=aa+bb/60
    if float(lon)<0:
        (c,d)=divmod(abs(float(lon)),100)
        cc=int(c)
        dd=float(d)
        lon_value=cc+(dd/60)
        lon_value=-lon_value
    else:
        (c,d)=divmod(float(lon),100)
        cc=int(c)
        dd=float(d)
        lon_value=cc+(dd/60)
    return lat_value, -lon_value   
  final_lat,final_lon=[],[]      
  for k in range(len(la)):
     [a,b]=dm2dd_single(la[k],lo[k])
     final_lat.append(a)
     final_lon.append(b)     
          
  return final_lat,final_lon
def whichArea(arg, lst):
    #Calculate certain point belongs to which area.
    i = len(lst)//2
    if i != 0: 
        if arg >= lst[i]:
            r = i + whichArea(arg, lst[i:])
        elif arg < lst[i]:
            r = whichArea(arg, lst[:i])
    else: r = i
    return r  
def draw_basemap(fig, ax, lonsize, latsize, interval_lon=2, interval_lat=2):
    ax = fig.sca(ax)
    dmap = Basemap(projection='cyl',
                   llcrnrlat=min(latsize)-0.01,
                   urcrnrlat=max(latsize)+0.01,
                   llcrnrlon=min(lonsize)-0.01,
                   urcrnrlon=max(lonsize)+0.01,
                   resolution='h',ax=ax)
    dmap.drawparallels(np.arange(int(min(latsize)),
                                 int(max(latsize))+1,interval_lat),
                       labels=[1,0,0,0], linewidth=0,fontsize=10)
    dmap.drawmeridians(np.arange(int(min(lonsize))-1,
                                 int(max(lonsize))+1,interval_lon),
                       labels=[0,0,0,1], linewidth=0,fontsize=10)
    dmap.drawcoastlines()
    dmap.fillcontinents(color='grey')
    dmap.drawmapboundary()
def stick_plot(time, u, v, **kw):
    width = kw.pop('width', 0.002)
    headwidth = kw.pop('headwidth', 0)
    headlength = kw.pop('headlength', 0)
    headaxislength = kw.pop('headaxislength', 0)
    angles = kw.pop('angles', 'uv')
    ax = kw.pop('ax', None)
    if angles != 'uv':
        raise AssertionError("Stickplot angles must be 'uv' so that"
              "if *U*==*V* the angle of the arrow on"
              "the plot is 45 degrees CCW from the *x*-axis.")
    time, u, v = map(np.asanyarray, (time, u, v))
    if not ax:
        fig, ax = plt.subplots()
    q = ax.quiver(date2num(time), [[0]*len(time)], u, v,angles='uv',
                  width=width,color='blue',**kw)
    ax.axes.get_yaxis().set_visible(False)
    dates = mpl.dates.drange(min(time),max(time),datetime.timedelta(days=1))
    dateFmt = mpl.dates.DateFormatter('%m/%d')
    ax.set_xticks(dates)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.xaxis.set_major_formatter(dateFmt)
    ax.xaxis.set_major_locator(mdates.DayLocator())
    return q
def getemolt_latlon(site):
    """
    get lat, lon, and depth for a particular emolt site 
    """
    import numpy as np
    #urllatlon = 'http://gisweb.wh.whoi.edu:8080/dods/whoi/emolt_site?emolt_site.SITE,emolt_site.LAT_DDMM,emolt_site.LON_DDMM,emolt_site.ORIGINAL_NAME,emolt_site.BTM_DEPTH&emolt_site.SITE='
    urllatlon = 'http://comet.nefsc.noaa.gov:8080/erddap/tabledap/eMOLT.csv?latitude,longitude,depth&SITE="'+str(site)+'"&distinct()'
    df=pd.read_csv(urllatlon,skiprows=[1])
    #dd=mean(df["depth"])
    dd=max(df["depth"])
    return df.latitude[0], df.longitude[0], dd
def str2list(s, bracket=False):
    '''
    a is a string converted from a list
    a = '[3,5,6,7,8]'
    b = str2list(a, bracket=True)
    or
    a = '3,4,5,6'
    b = str2list(a)
    '''
    if bracket:
        s = s[1:-1]
    s = s.split(',')
    s = [float(i) for i in s]
    return s
def str2ndlist(arg, bracket=False):
    '''
    convert list full of str to multidimensional arrays
    '''
    ret = []
    for i in arg:
        a = str2list(i, bracket=bracket)
        ret.append(a)
    # ret = np.array(ret)
    return ret

def getobs_tempsalt(site):
    """
    Function written by Jim Manning and used in "modvsobs"   
    get data from url, return datetime, temperature and depth
    """
    url = 'http://comet.nefsc.noaa.gov:8080/erddap/tabledap/eMOLT.csv?time,depth,sea_water_temperature&SITE="'+str(site)+'"&orderBy("time")'
    df=pd.read_csv(url,skiprows=[1])
    time=[]  
    for k in range(len(df)):
        time.append(datetime.datetime.strptime(df.time[k],'%Y-%m-%dT%H:%M:%SZ'))
    data={'time':time,'depth':df.depth.values,'temperature':df.sea_water_temperature.values}
    Data=pd.DataFrame(data)
    Data=Data.dropna()
    Data.index=range(len(Data))
    return Data
def get_emolt_data(alon,alat):
    'get nearest emolt site ,distance;then get fvcom depth of given site '
    SITES=['WP01','RS01','NM01','AC02','SJ02','TA24','RG01','DMF0','EP01',
           'BT03','MJ09','ON04','AP01','OM03','IC01','ES01','JV02','JP22',
           'KM01','PM04','BL01','ON02','PB02','BT04','AC01','PB03','MO01',
           'PM03','JP19','DM03','SR08','CW01','JP02','OD08','CC01','RR01',
           'PM02','PC01','NL02','DS01','PM01','SG01','JD01','RP01','MW03',
           'MC02','TS02','BT01','TA14','DMF8','RM03','JP16','JA02','JA01',
           'OM01','ET02','CP02','MA10','DMF9','WL01','JM01','RB01','MW01',
           'WD02','PF01','SJ01','LC01','BBHR','TH01','TS01','PW01','RB02',
           'GR01','GS01','RA01','GR02','OC01','DK01','BD01','JT04','MF02',
           'CJ01','RM04','AB01','BS02','DC01','MM01','WD01','BI02','BI01',
           'JC01','BC01','CP01','KO01','BM02','RM02','BF01','RM01','ET01',
           'BA03','DMF6','NL01','DMF5','DJ01','BM01','JS02','DMF2','DMF7', 
           'JS06','DMF4','DMF1','WHAQ','BN01','TA15','BA02','AG01','BA01']    #number of these sites` observations are >10000.  
    url='http://comet.nefsc.noaa.gov:8080/erddap/tabledap/eMOLT.csv?SITE,latitude,longitude,depth&distinct()'
    df=pd.read_csv(url,skiprows=[1])
    sites=pd.Series(df.SITE)
    lons=pd.Series(df.longitude)
    lats=pd.Series(df.latitude)
    depths=pd.Series(df.depth)      #get data from website
    indx=[]    
    for i in SITES:
        for j in range(len(df)):
            if sites[j]==i:
                indx.append(j)     #get above sites` index
    data={'site':sites[indx].values,'lon':lons[indx].values,'lat':lats[indx].values}
    Data=pd.DataFrame(data,index=range(len(indx)))   #get above sites` data
    dist=DIST(alon,alat,Data['lon'].values,Data['lat'].values)   
    nearest_index=np.argmin(dist)        #calculate nearest index
    data_nearest=getobs_tempsalt(Data['site'][nearest_index])   #get nearest site data
    return dist[nearest_index],data_nearest,Data['site'][nearest_index]
def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr
def shrink(a,b):
    """Return array shrunk to fit a specified shape by triming or averaging.
    a = shrink(array, shape)
    array is an numpy ndarray, and shape is a tuple (e.g., from
    array.shape). a is the input array shrunk such that its maximum
    dimensions are given by shape. If shape has more dimensions than
    array, the last dimensions of shape are fit.
    as, bs = shrink(a, b)
    If the second argument is also an array, both a and b are shrunk to
    the dimensions of each other. The input arrays must have the same
    number of dimensions, and the resulting arrays will have the same
    shape.
    """
    if isinstance(b, np.ndarray):
        if not len(a.shape) == len(b.shape):
            raise(Exception, \
                  'input arrays must have the same number of dimensions')
        a = shrink(a,b.shape)
        b = shrink(b,a.shape)
        return (a, b)
    if isinstance(b, int):
        b = (b,)
    if len(a.shape) == 1:                # 1D array is a special case
        dim = b[-1]
        while a.shape[0] > dim:          # only shrink a
            if (dim - a.shape[0]) >= 2:  # trim off edges evenly
                a = a[1:-1]
            else:                        # or average adjacent cells
                a = 0.5*(a[1:] + a[:-1])
    else:
        for dim_idx in range(-(len(a.shape)),0):
            dim = b[dim_idx]
            a = a.swapaxes(0,dim_idx)        # put working dim first
            while a.shape[0] > dim:          # only shrink a
                if (a.shape[0] - dim) >= 2:  # trim off edges evenly
                    a = a[1:-1,:]
                if (a.shape[0] - dim) == 1:  # or average adjacent cells
                    a = 0.5*(a[1:,:] + a[:-1,:])
            a = a.swapaxes(0,dim_idx)        # swap working dim back
    return a
def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr
def bbox2ij(lon,lat,bbox):
    """Return indices for i,j that will completely cover the specified bounding box.     
    lon,lat = 2D arrays that are the target of the subset
    bbox = list containing the bounding box: [lon_min, lon_max, lat_min, lat_max]
    """
    bbox=np.array(bbox)
    mypath=np.array([bbox[[0,1,1,0]],bbox[[2,2,3,3]]]).T
    p = path.Path(mypath)
    points = np.vstack((lon.flatten(),lat.flatten())).T   
    n,m = np.shape(lon)
    inside=[]
    for i in range(len(points)):
        inside.append(p.contains_point(points[i]))
    inside = np.array(inside).reshape((n,m))
    ii,jj = np.meshgrid(range(m),range(n))
    inside_i,inside_j=[],[]
    for i in range(len(inside)):
        for j in range(len(inside[i])):
            if inside[i][j]==1:
                inside_i.append(ii[i][j])
                inside_j.append(jj[i][j])
    return min(inside_i),max(inside_i),min(inside_j),max(inside_j)
def gmt_to_loc(utc):
    from_zone = tz.gettz('UTC')
    to_zone = tz.gettz('America/New_York')
    from_zone = tz.tzutc()
    to_zone = tz.tzlocal()
    utc = utc.replace(tzinfo=from_zone)
    central = utc.astimezone(to_zone)
    return central
def get_roms_url(time,method):
    '''
    use time to get roms url.     
    time is datetime.
    method is forecast or hindcast. 
    if method=='forecast':
        time is datetime
    if method=='hindcast':
        time is [datetime0,datetime1]
    '''
    if method=='forecast':
        date=time-dt.timedelta(days=3)
        str_date=dt.datetime.strftime(date.date(),"%Y-%m-%d")
        url='http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/his/runs/ESPRESSO_Real-Time_v2_History_RUN_'+str_date+'T00:00:00Z'
    if method=='hindcast':
        start=time[0]
        end=time[1]
        url=[]
        diff_days=(end.date()-start.date()).days
        for i in range(diff_days):
            str_date=dt.datetime.strftime((start+timedelta(days=i)).date(),"%Y%m%d")
            url_each='http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/his/files/espresso_his_'+str_date+'_0000_0001.nc'
            url.append(url_each)
    return url    
