#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 16:52:23 2018

@author: jordianguera
"""

import h5py
import os
import numpy as np
import matplotlib.pyplot as plt

ndays = 360
stations = ['FR.CALF.00.HHZ', 'FR.EILF.00.HHZ', 'FR.ESCA.01.HHZ', 'FR.MON.00.HHZ', 'FR.MVIF.00.HHZ', 'FR.PRIMA.00.HHZ', 'FR.SAOF.00.HHZ']
rmax = 7
init = 2016183
seven_days = []
seven_working_stations = []
f = []

for i in range(ndays):
    
    if (init + i <= 2016366):
        file_h5 = str(init + i)
    else:
        file_h5 = str(2016817 + i)
    
    dt = []
    available = []

    os.chdir('/Users/jordianguera/Desktop/ECMI_Data/dt/AlpArray')
    filename_h5 = file_h5 + '.h5'
    f = h5py.File(filename_h5, 'r')    
    file = list(f.keys())[0]
    data = list(f[file])    
    
    os.chdir('/Users/jordianguera/Desktop/ECMI_Data/dt/AlpArray/METRICS')  
    filename_stations = file_h5 + "_local_pair_dist.txt"
    
    with open(filename_stations) as x:
        m = 0
        for line in x:            
            g = line.split(" ")
            for k in range(rmax):
                if (g[0] == stations[k]):
                    for l in range(rmax):
                        if (g[1] == stations[l]):
                            available.append(g)
                            
                            dt.append(data[m])
            m += 1      

    seven_days.append(dt)          
    seven_working_stations.append(available)
            
#%% ACTIVE STATIONS

seven_active = [None] * ndays
for i in range(ndays):
    seven_active[i] = [None] * len(seven_working_stations[i])
    for j in range(len(seven_working_stations[i])):
        seven_active[i][j] = [None] * 2
        seven_active[i][j][0] = seven_working_stations[i][j][0]
        seven_active[i][j][1] = seven_working_stations[i][j][1]
        
#%% BUILD CONNECTIONS ARRAY

seven_stations_connections = []

for i in range(rmax-1):
    for j in range(i+1, rmax):
        name = [stations[i], stations[j]]
        seven_stations_connections.append(name)

#%% BUILD MATRIX        

#Initialise dt_matrix  
seven_dt_matrix = [None] * len(seven_stations_connections)
for i in range(len(seven_stations_connections)):
    seven_dt_matrix[i] = [0] * (24*ndays)
 
#Fill dt_matrix
for day in range(ndays):
    for i in range(len(seven_stations_connections)):
        #print(day, i)
        for station in range(len(seven_active[day])):
            if((seven_stations_connections[i][0] == seven_active[day][station][0]) & (seven_stations_connections[i][1] == seven_active[day][station][1])):
                for hour in range(day*24, (day+1)*24):
                    #print(station, hour)
                    seven_dt_matrix[i][hour] = seven_days[day][station][0][hour - day*24]
                    

 
        
        
        
        
        
        
        
        
