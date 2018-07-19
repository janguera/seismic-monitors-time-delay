import h5py
import os
import matplotlib.pyplot as plt

ndays = 360
stations = ['FR.CALF.00.HHZ', 'FR.EILF.00.HHZ', 'FR.ESCA.01.HHZ', 'FR.MON.00.HHZ', 'FR.MVIF.00.HHZ', 'FR.PRIMA.00.HHZ', 'FR.SAOF.00.HHZ']
all_stations = ['FR.ARBF.00.HHZ', 	'FR.ARTF.00.HHZ', 	'FR.BANN.00.HHZ', 	'FR.BLAF.00.HHZ', 	'FR.BSTF.00.HHZ', 	'FR.CALF.00.HHZ', 	'FR.CFF.00.HHZ', 	'FR.CLAF.00.HHZ', 	'FR.COLF.00.HHZ', 	'FR.EILF.00.HHZ', 	'FR.ENAUX.00.HHZ', 	'FR.ESCA.00.HHZ', 	'FR.FLAF.00.HHZ', 	'FR.GRN.00.HHZ', 	'FR.ISO.01.HHZ', 	'FR.LBL.01.HHZ', 	'FR.MLYF.00.HHZ', 	'FR.MON.00.HHZ', 	'FR.MVIF.00.HHZ', 	'FR.OGAG.00.HHZ', 	'FR.OGCB.00.HHZ', 	'FR.OGDI.00.HHZ', 	'FR.OGGM.00.HHZ', 	'FR.OGMO.00.HHZ', 	'FR.OGMY.00.HHZ', 	'FR.OGS1.00.HHZ', 	'FR.OGS2.00.HHZ', 	'FR.OGS3.00.HHZ', 	'FR.OGSA.00.HHZ', 	'FR.OGSM.00.HHZ', 	'FR.OGVG.00.HHZ', 	'FR.PLDF.00.HHZ', 	'FR.PRIMA.00.HHZ', 	'FR.RSL.00.HHZ', 	'FR.RUSF.00.HHZ', 	'FR.SAOF.00.HHZ', 	'FR.SAUF.00.HHZ', 	'FR.SPIF.00.HHZ', 	'FR.SURF.00.HHZ', 	'FR.TRBF.00.HHZ', 	'FR.TRIGF.00.HHZ', 	'FR.TURF.00.HHZ', 	'FR.VAL4.06.HHZ', 	'FR.SSB.05.HHZ', 	'FR.A174A.03.HHZ', 	'FR.A175A.07.HHZ', 	'FR.A176A.01.HHZ', 	'FR.A177A.04.HHZ', 	'FR.A178A.05.HHZ', 	'FR.A180A.00.HHZ', 	'FR.A181A.00.HHZ', 	'FR.A183A.00.HHZ', 	'FR.A184A.00.HHZ', 	'FR.A185A.00.HHZ', 	'FR.A186A.00.HHZ', 	'FR.A186B.00.HHZ', 	'FR.A187A.00.HHZ', 	'G.A188A.10.HHZ', 	'Z3.A191A.00.HHZ', 	'Z3.A192A.00.HHZ', 	'Z3.A193A.00.HHZ', 	'Z3.A194A.00.HHZ', 	'Z3.A196A.00.HHZ', 	'Z3.A199A.00.HHZ', 	'Z3.A200A.00.HHZ', 	'Z3.A201A.00.HHZ', 	'Z3.A202A.00.HHZ', 	'Z3.A204A.00.HHZ', 	'Z3.A205A.00.HHZ', 	'Z3.A206A.00.HHZ', 	'Z3.A215A.00.HHZ', 	'Z3.A216A.00.HHZ', 	'Z3.A217A.00.HHZ', ]
rmax = 73
init = 2016183
days = []
working_stations = []
f = []
j = 0

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
                if (g[0] == all_stations[k]):
                    for l in range(rmax):
                        if (g[1] == all_stations[l]):
                            available.append(g)
                            
                            dt.append(data[m])

                            j += 1
            m += 1      

    days.append(dt)          
    working_stations.append(available)

#%% BUILD CONNECTIONS ARRAY

stations_connections = []

for i in range(rmax-1):
    for j in range(i+1, rmax):
        name = [all_stations[i], all_stations[j]]
        stations_connections.append(name)


#%% ACTIVE STATIONS

active = [None] * ndays
for i in range(ndays):
    active[i] = [None] * len(working_stations[i])
    for j in range(len(working_stations[i])):
        active[i][j] = [None] * 2
        active[i][j][0] = working_stations[i][j][0]
        active[i][j][1] = working_stations[i][j][1]
        
#%% BUILD MATRIX        

#Initialise dt_matrix  
dt_matrix = [None] * len(stations_connections)
for i in range(len(stations_connections)):
    dt_matrix[i] = [0] * (24*ndays)
    
#Fill dt_matrix
for day in range(ndays):
    print(day)
    for i in range(len(stations_connections)):
        #print(day, i)
        for station in range(len(active[day])):
            if((stations_connections[i][0] == active[day][station][0]) & (stations_connections[i][1] == active[day][station][1])):
                for hour in range(day*24, (day+1)*24):
                    #print(station, hour)
                    dt_matrix[i][hour] = days[day][station][0][hour - day*24]
   
#%% SAVE THE dt MATRIX

from scipy import sparse
import numpy

a = numpy.asarray(dt_matrix)
CSC_dt = sparse.csc_matrix(a)
print(CSC_dt)
sparse.save_npz('CSC_dt.npz', CSC_dt)
                 
#%% NUMBER OF WORKING CONNECTIONS

nconnect = []
for i in range(ndays):
    nconnect.append(len(working_stations[i]))   
    #print(len(working_stations[i]))
plt.plot(nconnect)    
plt.title("Number of working connections")
plt.show()


#%% WORKING STATIONS

n_working = [None] * ndays
for day in range(ndays):
    stations_working = [0] * rmax 
    for j in range(len(active[day])):
        for k in range(rmax):
            if((all_stations[k] == active[day][j][0]) | (all_stations[k] == active[day][j][1])):
                stations_working[k] = 1
    n_working[day] = sum(stations_working)



plt.plot(n_working)    
plt.title("Number of working seismic monitors")
plt.show()

#%% WORKING STATIONS-CONNECTIONS CORRELATION

nconnect_std = []
nworking_std = []

for day in range(ndays):
    nconnect_std.append(nconnect[day]/max(nconnect))
    nworking_std.append(n_working[day]/max(n_working))
    
plt.plot(nconnect_std)    
plt.title("Standarised number of working connections")
plt.show()

plt.plot(nworking_std)    
plt.title("Standarised number of working seismic monitors")
plt.show()

correlation_stations = []

for day in range(ndays):
    correlation_stations.append(nconnect[day]/n_working[day])

plt.plot(correlation_stations)    
plt.title("Number of connections per working seismic monitor")
plt.show()

#%% DAILY dt sum

daily_sum = []
N_connections = len(stations_connections)
for hour in range(ndays*24):
    sum = 0
    for i in range(N_connections):
        sum += dt_matrix[i][hour]
    daily_sum.append(sum)
    
daily_avg = []    
for hour in range(ndays*24):
    daily_avg.append(daily_sum[hour]/N_connections)
    
plt.plot(daily_avg)    
plt.title("Daily dt average")
plt.show()
    
plt.plot(daily_sum)    
plt.title("Daily dt sum")
plt.show()
