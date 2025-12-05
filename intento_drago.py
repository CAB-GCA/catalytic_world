from fun_gilles import *
import matplotlib.pyplot as plt
import numpy as np
import pickle

file = 'examples/reactions_autocat.txt'
reactions = read_file(file)
species = obtain_species(reactions)

n_iterations = 5e5
method = "Protocell" # Gillespie or Deterministic
# Reaction constants:
k = [1]*8 # len(k)= # de reacciones
# Volume:
V = 1000
initial_food = [1000]*4 + [500]*2 + [0]*2 # initial molecules number
final_volume = []
k_change = np.logspace(-3,4,8)

for k_i in k_change:
    k[2] = k_i
    print(f"Performing simulation for k = {k_i}")
    abundances, times, volumes = chemistry(method, n_iterations, reactions,
                            initial_food, k, V)
    print(f"Simulation ended. Final volume = {volumes[-1]}")
    
    try:
        print(f"Volume at t= 10 = {volumes[times>10][0]}")
        final_volume.append(volumes[times>10][0])
    except IndexError:
        final_volume.append(volumes[-1])
        
    data = {
        'abundance':abundances,
        'times':times,
        'volumes':volumes
    }
    
    with open(f'data_for_{k_i}.txt', 'wb') as data_file:
        pickle.dump(data, data_file)
    
    