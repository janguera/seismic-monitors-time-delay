import h5py
import os
import matplotlib.pyplot as plt
import numpy as np

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
    file = list(f.keys())[2]
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
    for i in range(len(stations_connections)):
        #print(day, i)
        for station in range(len(active[day])):
            if(((stations_connections[i][0] == active[day][station][0]) & (stations_connections[i][1] == active[day][station][1]) | (stations_connections[i][1] == active[day][station][0]) & (stations_connections[i][0] == active[day][station][1]))):
                for hour in range(day*24, (day+1)*24):
                    #print(station, hour)
                    dt_matrix[i][hour] = days[day][station][0][hour - day*24]
   
#%% SAVE THE dt MATRIX

from scipy import sparse

dt_np_matrix = np.asarray(dt_matrix)
CSC_dt = sparse.csc_matrix(dt_np_matrix)
sparse.save_npz('CSC_dt.npz', CSC_dt)
                 
#%% NUMBER OF WORKING CONNECTIONS

nconnect = []
for i in range(ndays):
    nconnect.append(len(working_stations[i]))   
    #print(len(working_stations[i]))
    
plt.plot(nconnect)    
plt.title("Number of working connections", fontsize=14)
plt.xlabel('Day', fontsize=12)
plt.ylabel('Connections', fontsize=12)
plt.show()

#%% NUMBER OF WORKING STATIONS

n_working = []
for day in range(ndays):
    stations_working = [0] * rmax 
    for j in range(len(active[day])):
        for k in range(rmax):
            if((all_stations[k] == active[day][j][0]) | (all_stations[k] == active[day][j][1])):
                stations_working[k] = 1
    n_working.append(np.sum(stations_working))



plt.plot(n_working)    
plt.title("Number of working seismic monitors", fontsize=14)
plt.xlabel('Day', fontsize=12)
plt.ylabel('Active monitors', fontsize=12)
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
plt.title("Average daily connections per working seismic monitor", fontsize=14)
plt.xlabel('Day', fontsize=12)
plt.ylabel('Connections per monitor', fontsize=12)
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
    daily_avg.append(daily_sum[hour]/active)

plt.figure(figsize=(20,10))    
plt.plot(daily_avg)    
plt.title("Daily dt average")
plt.show()
    
plt.plot(daily_sum)    
plt.title("Daily dt sum")
plt.show()

#%% SINGULAR VALUE DECOMPOSITION

#Here we can see that the relative importance of connections within a day is more 
#or less the same within the active ones

#Day 1
day = 0
day1_matrix = [None] * len(working_stations[day])

for i in range(len(working_stations[day])):
    day1_matrix[i] = [None] * 24
    for j in range(24):
        day1_matrix[i][j] = dt_matrix[i][j]

s1 = np.linalg.svd(day1_matrix, full_matrices=True, compute_uv=False)

#Day 179
day = 180
day2_matrix = [None] * len(working_stations[day])

for i in range(len(working_stations[day])):
    day2_matrix[i] = [None] * 24
    for j in range(24):
        day2_matrix[i][j] = dt_matrix[i][j+day]

s2 = np.linalg.svd(day2_matrix, full_matrices=True, compute_uv=False)

#Day 280
day = 279
day3_matrix = [None] * len(working_stations[day])

for i in range(len(working_stations[day])):
    day3_matrix[i] = [None] * 24
    for j in range(24):
        day3_matrix[i][j] = dt_matrix[i][j+day]

s3 = np.linalg.svd(day3_matrix, full_matrices=True, compute_uv=False)

plt.plot(s1)
plt.plot(s2)
plt.plot(s3)
plt.show()


#Binary full matrix - we will identify the index of the eigenvector which has the 
#most relative importance (only by number of connections) - we do this because we 
#couln't come up with a conclusion on a daily basis, so we take the whole data

binary_matrix = [0] * len(stations_connections)

for i in range(len(stations_connections)):
    binary_matrix[i] = [0] * (24*ndays)
    for j in range(24*ndays):
        if(dt_matrix[i][j] > 0):
            binary_matrix[i][j] = 1
            
u_bin, s_bin, v_bin = np.linalg.svd(binary_matrix, full_matrices=True, compute_uv=True)

plt.title("Singular values for the binary connection matrix", fontsize=14)
plt.xlabel('Eigenvector\'s index', fontsize=12)
plt.ylabel('Relative weight', fontsize=12)
plt.plot(s_bin)
plt.show()

plt.title("Singular values for the binary connection matrix", fontsize=14)
plt.xlabel('Eigenvector\'s index', fontsize=12)
plt.ylabel('Relative weight', fontsize=12)
plt.plot(v_bin[0])
plt.show()

#Once we know the index for which the eigenvector is the most important we will plot 
#the relative weight by number of connections to see what connections have the largest
#relative weight and what is the limit at which the connections don't have any more influence 


u, s, v = np.linalg.svd(dt_matrix, full_matrices=True, compute_uv=True)

plt.title("Graph's connections influence",fontsize=16)
plt.xlabel('Number of connections', fontsize=12)
plt.ylabel('Relative weight', fontsize=12)
plt.plot(u[0])
plt.show()

#%% FINDING THE MOST IMPORTANT STATIONS

#Most influent connections

u_np = np.array(u[0])

sort_index = np.argsort(u_np)
u_sorted = []
for i in range(len(u_np)):
    u_sorted.append(u_np[sort_index[i]])

plt.plot(u_sorted)

#We set a certain threshold to detect the most influentiable connections
threshold = 0.05
important_connections = []
for i in range(len(u_np)):
    if ((u_sorted[i] <= -threshold) | (u_sorted[i] >= threshold)):
        important_connections.append(sort_index[i])
        
imp_conn_names = []
for i in range(len(important_connections)):
    imp_conn_names.append(stations_connections[important_connections[i]][0])
    imp_conn_names.append(stations_connections[important_connections[i]][1])

#Here we identify the stations with the most influentiable connections         
unique_stations = list(set(imp_conn_names))

count = [0] * len(unique_stations)
for i in range(len(unique_stations)):
    for j in range(len(imp_conn_names)):
        if (unique_stations[i] == imp_conn_names[j]):
            count[i] += 1
print(count[42])        
print(unique_stations[6])        
print(unique_stations[20])        
print(unique_stations[42])        

sort_unique_stations = list(np.sort(count))
print(sort_unique_stations[0])    
    