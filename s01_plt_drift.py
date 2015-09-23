# -*- coding: utf-8 -*-
"""
Created on Mon Apr 09 23:05:10 2012

@author: Vitalii Sheremet, FATE Project
"""
# NEFSC Drifter Archive analysis
# ipython --pylab run s01_plt_drift.py
# modified 2015-09-23 new dowload and processing approach

#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime


#/home/vsheremet/FATE/DriftFVCOM2/drifters_archive.csv
#id,time,latitude,longitude,depth,sea_water_temperature
#,UTC,degrees_north,degrees_east,m,degree_C
#100230811,2010-12-17T17:00:00Z,32.3293,-77.9628,-1.0,NaN
#45387,2004-06-29T12:00:00Z,43.6422,-68.9547,15.0,10.13
FN='drifters_archive.csv' # whole archive ~30MB in 2015
D = np.genfromtxt(FN,dtype=None,names=['id','time','latitude','longitude','depth','sea_water_temperature'],delimiter=',',skip_header=2)

print('All fixes:'+str(len(D['id'])))
IDs = list(set(D['id']))
print('All tracks, unique drifter IDs:' + str(len(IDs)))

#t0=datetime.strptime('1990-01-01T00:00:00Z','%Y-%m-%dT%H:%M:%S.%f')
t0=datetime.strptime('1990-01-01T00:00:00Z','%Y-%m-%dT%H:%M:%SZ')

k=0    
while k in range(len(IDs)):
# browse through ids     
    ID=IDs[k]
    print('ID:',ID)
    i=np.argwhere(D['id']==ID).flatten()

# single drifter track
    T=D['time'][i]  # string
    lon=D['longitude'][i]
    lat=D['latitude'][i]
    t = np.array([datetime.strptime(T[kt],'%Y-%m-%dT%H:%M:%SZ') for kt in range(len(T))]) # datetime
#    tp = np.array([(t[kt]-t0).total_seconds()/60/60/24. for kt in range(len(t))]) 

    plt.figure(1)    
    plt.title(str(ID))
    plt.subplot(211)
    plt.plot(t,lat,'b-')
    plt.xlabel('time, UTC')
    plt.ylabel('lat')
    plt.title(str(ID))
    plt.grid(True)

    plt.subplot(212)
    plt.plot(t,lon,'r-')
    plt.xlabel('time, UTC')
    plt.ylabel('lon')
    plt.grid(True)
    plt.show()
    
    plt.figure(2)
    #trdh=(trd[1:]+trd[:-1])*0.5
    lath=(lat[1:]+lon[:-1])*0.5
    clat=np.cos(lath*np.pi/180.)
    d=np.diff(t)
    dt=np.array([d[kt].total_seconds() for kt in range(len(d))]).flatten() 
   
    vlon=np.diff(lon)*111111.*clat/dt
    vlat=np.diff(lat)*111111./dt
    plt.plot(t[:-1],vlon,'r.-',t[:-1],vlat,'b.-')
    plt.xlabel('time UTC')
    plt.ylabel('vlon (r), vlat (b) m/s')
    plt.title(ID)
    plt.grid(True)
    plt.show()
    
    plt.figure(3)
    plt.plot(lon,lat,'r.-',lon[0],lat[0],'go',lon[-1],lat[-1],'b*')
    plt.xlabel('lon')
    plt.ylabel('lat')
    plt.title(ID)
    plt.grid(True)
    plt.show()
    

    print(str(k)+'  '+str(ID))
    print('Enter: b - step 1 back, # [0:'+str(len(IDs)-1)+'] - go to #, anything else - step 1 forward')    
    a=raw_input()
    
    if a=='b':
        k=k-1  # step back
    else:
        try:
            k=int(a)
        except:
            k=k+1  # step forward by default
            
    plt.close(1)
    plt.close(2)
    plt.close(3)

