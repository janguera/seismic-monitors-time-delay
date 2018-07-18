import numpy as np
def graph_denoising(noisy, measure_estimates, epsilon, learning_rate, max_iterations):
    t = 0
    W = noisy
    E = np.sum(np.power((graph_measures(W) - measure_estimates), 2))
    while E > epsilon and t < max_iterations:
        measures = graph_measures(W)
        measure_derivatives = graph_measure_derivatives(W)
        measures_minus_estimates_times_derivatives = [
            np.dot((measures[i] - measure_estimates[i]), measure_derivatives) for i in range(len(measure_estimates))
        ]
        W = W - learning_rate * np.sum(measures_minus_estimates_times_derivatives, axis=0)
        W[W < 0] = 0
        W[W > 1] = 1
        E = np.sum(np.power((graph_measures(W) - measure_estimates), 2))
        t = t + 1
        if t >= max_iterations:
            print("The algorithm stopped because of reaching max iterations.")
    return W, t

def graph_measures(W):
    # Fill the graph measures here
    # return a list of length n, where n is the number of measures we are using
    return np.array([1,2,3])

def graph_measure_derivatives(W):
    # Fill the graph measure derivatives here
    # return something with the shape of (W,n) where n is the number of measures we are using
    return np.array([W, W, W])

a = np.array([[0,0.5],[0.4,0.6]])
estimates = [1,2,3]
epsilon = 100
learning_rate = 0.05

graph_denoising(a, estimates, epsilon, learning_rate)