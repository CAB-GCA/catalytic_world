import numpy as np
import sys
sys.path.append('..')  # allow importing from parent directory
from fun_gilles import *
import pickle
import pandas as pd

f = 'examples/reactions_autocat.txt'
reactions = read_file(f)
species = obtain_species(reactions)

# --- BARRIDO DE CONDICIONES INICIALES ---
n_iterations = 5e5
method = "Protocell" # Gillespie or Deterministic
# Reaction constants:
k = [1]*8 # len(k)= # de reacciones
# Volume:
V = 1000
initial_food = [1000]*4 + [0]*4 # initial molecules number


# cat_abundance = np.concatenate((np.linspace(0,1000,9, endpoint=True), np.linspace(1500, 7000, 12)))
# cat_abundance[0] = 1
# print(cat_abundance)

# simulations={}

# for initial_cat in cat_abundance:
#     initial_food[4] = initial_cat
#     abundances, times, volumes = chemistry(method, n_iterations, f,
#                                 initial_food, k, V)
#     print(f"Simulation completed for [ab]0 = {initial_cat/V}")
#     simulations[initial_cat/V]=(abundances, times, volumes)
    
# with open("barrido_ab.txt", 'ab') as file:
#     pickle.dump(simulations, file)
    
# # --- BARRIDO DE K_0 ---

# # Reaction constants:
# k = [1]*8 # len(k)= # de reacciones
# # Volume:
# V = 1000
# initial_food = [1000]*4 + [500]*2 + [0]*2 # initial molecules number
# k_change = np.logspace(-3,4,15)
# simulations_k0 = {}

# for k_i in k_change:
#     k[0] = k_i
#     abundances, times, volumes = chemistry(method, n_iterations, f,
#                             initial_food, k, V)
#     print(f"Simulation completed for k_0 = {k_i}")
#     simulations_k0[k_i] = (abundances, times, volumes)

# with open("barrido_k0.txt", "ab") as file:
#     pickle.dump(simulations_k0, file)

# --- BARRIDO DE K_2 ---

# Reaction constants:
k = [1]*8 # len(k)= # de reacciones
# Volume:
V = 1000
initial_food = [1000]*4 + [500]*2 + [0]*2 # initial molecules number
final_volume = []
k_change = np.logspace(-3,4,15)
simulations_k2 = {}

for k_i in k_change:
    k[2] = k_i
    abundances, times, volumes = chemistry(method, n_iterations, f,
                            initial_food, k, V)
    print(f"Simulation completed for k_2 = {k_i}")

    simulations_k2[k_i] = (abundances, times, volumes)
    
with open("barrido_k2.txt", "ab") as file:
    pickle.dump(simulations_k2, file)