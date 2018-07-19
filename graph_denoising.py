import numpy as np
import pandas as pd
np.set_printoptions(precision=2)

def graph_denoising(noisy, measure_estimates, epsilon, learning_rate, max_iterations):
    t = 0
    # print(np.shape(noisy))
    W = noisy
    # print(graph_measures(W))
    E = np.sum(np.power(np.subtract(graph_measures(W), measure_estimates), 2))
    # print(graph_measures(W) - measure_estimates)
    #print(E)
    while E > epsilon and t < max_iterations:
        measures = graph_measures(W)
        # print("shape of measures " + str(np.shape(measures)))
        measure_derivatives = graph_measure_derivatives(W)
        # print("shape of derivatives " +str(np.shape(measure_derivatives)))
        #print(len(measure_estimates))
        #print(len(measures))
        # print(measures[0])
        # print(measure_derivatives[0])
        measures_minus_estimates_times_derivatives = [
            np.dot((measures[i] - measure_estimates[i]), measure_derivatives[i]) for i in range(len(measure_estimates))
        ]
        # print(measures_minus_estimates_times_derivatives)
        # print(np.shape(measures_minus_estimates_times_derivatives))
        # print(np.shape(np.sum(measures_minus_estimates_times_derivatives, axis=0)))
        W = W - learning_rate * np.sum(measures_minus_estimates_times_derivatives, axis=0)
        # print(np.shape(W))
        W[W < 0] = 0
        W[W > 1] = 1
        E = np.sum(np.power((graph_measures(W) - measure_estimates), 2))
        # print(E)
        t = t + 1
        if t >= max_iterations:
            print("The algorithm stopped because of reaching max. Error at the end was " + str(E))
    return W, t

def graph_measures(W):
    # Fill the graph measures here
    # return a list of length n, where n is the number of measures we are using
    measures = []
    for i in range(len(W)):
        measures.append(distance_measure(W, i))

    return np.array(measures)

def distance_measure(W, i):
    return np.sum(W, axis=1)[i] # = [1,2,3,4,5,6,7][i]


def graph_measure_derivatives(W):
    # Fill the graph measure derivatives here
    # return something with the shape of (W,n) where n is the number of measures we are using
    return np.array([np.ones(np.shape(W)) for i in range(len(W))])

##### determinig the W_d.
filename_info = 'C:/MyTemp/ECMI_Data/dt/AlpArray/METRICS/2016366_local_pair_dist.txt'

distances = np.zeros((7,7))
# print(distances)
file = open(filename_info)
names = ["FR.CALF.00.HHZ",
         "FR.EILF.00.HHZ",
         "FR.ESCA.01.HHZ",
         "FR.MON.00.HHZ",
         "FR.MVIF.00.HHZ",
         "FR.PRIMA.00.HHZ",
         "FR.SAOF.00.HHZ"]
i = 0
for line in file:
    line_split = line.split(" ")
    if line_split[0] in names and line_split[1] in names:
        distances[names.index(line_split[0]), names.index(line_split[1])] = line_split[2]
        distances[names.index(line_split[1]), names.index(line_split[0])] = line_split[2]
        i = i+1
file.close()
# print(distances)
for i in range(7):
    distances[i,i] = np.nan
one_by_r_distances = 1/distances
# print(one_by_r_distances)
max_one_by_r = np.nanmax(one_by_r_distances)
min_one_by_r = np.nanmin(one_by_r_distances)
# print(max_one_by_r)
normalized_one = one_by_r_distances /max_one_by_r
normalized_one = np.nan_to_num(normalized_one, 0)
normalized_zero_to_one = (one_by_r_distances - min_one_by_r) / (max_one_by_r - min_one_by_r)
normalized_zero_to_one = np.nan_to_num(normalized_zero_to_one, 0)
W_d = normalized_zero_to_one
#####

#### read W_t
array = pd.read_csv("C:/MyTemp/ECMI_Data/dt/AlpArray/day_1/dt_day1.csv", delimiter=";").values[:,1].astype(np.float32)
# array = np.array
result = np.zeros((7,7))
#result[0, 1:7] = array[0:6]
#result[1, 0:1] = array[6]
#result[1, 2:7] = array[7:11]
#result[2, 0:2] = array[11:13]
#result[2, 3:7] = array[13:16]

i = 0
j = 0
for t in range(len(array)):
    if j == i:
        j = j + 1
    #print("t= " + str(t))
    #print(array[t])
    result[i, j] = array[t]
    result[j, i] = array[t]
    # print(t)
    #print("i = " + str(i))
    #print("j = " + str(j))
    if j < 6:
        j = j + 1
    else:
        j = i + 1
        i = i + 1
W_t = result

####

# W_t = np.ones((7,7)) # time delay matrix
# W_d = np.array([[0,0.4],[0.8,0.2]]) # distance matrix
estimates = graph_measures(W_d) #estimates
epsilon = 1 # This should be quite small
learning_rate = 0.05 # I have no idea what this should be, we have to try different values

W, t = graph_denoising(W_t, estimates, epsilon, learning_rate, max_iterations=1000000) # ready to run when we get the matrices.
print(W_t)
print(W)
