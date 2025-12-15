import numpy as np
import sys
sys.path.append('..')  # allow importing from parent directory
from fun_gilles import *
import pickle
import pandas as pd

def barrido_k(k_n, k_values):
    # Reaction constants:
    k = [1]*8 # len(k)= # de reacciones
    # Volume:
    V = 1000
    initial_food = [1000]*4 + [500]*2 + [0]*2 # initial molecules number
    final_volume = []
    output_file = f"barrido_k{k_n}_stream.pkl"

    try: 
        with open(output_file, "ab") as file:

            for k_i in k_values:
                k[k_n] = k_i
                print(f"Starting simulation for k_2 = {k_i}")
                abundances, times, volumes = chemistry(method, n_iterations, f,
                                        initial_food, k, V)
                print(f"Simulation completed for k_2 = {k_i}")

                results_k2 = {k_i:(abundances, times, volumes)}
                pickle.dump(results_k2, file)

                file.flush()
    except Exception as e:
        print(f"An error has occured: {e}")

    finally: 
        print("Simulation run completed")

f = 'examples/reactions_autocat.txt'
reactions = read_file(f)
species = obtain_species(reactions)

# --- BARRIDO DE CONDICIONES INICIALES ---
n_iterations = 2e5
method = "Protocell" # Gillespie or Deterministic
# Reaction constants:
k = [1]*8 # len(k)= # de reacciones
k[0] = 1e2
# Volume:
V = 1000
initial_food = [1000]*4 + [0]*4 # initial molecules number


cat_abundance = (np.linspace(7000, 27000, 41))
print(cat_abundance)
output_file = "barrido_ab0_k0_e2.pkl"

try:
    with open(output_file, "ab") as file:

        for initial_cat in cat_abundance:
            initial_food[4] = initial_cat
            abundances, times, volumes = chemistry(method, n_iterations, f,
                                        initial_food, k, V)
            print(f"Simulation completed for [ab]0 = {initial_cat/V}")
            results = {initial_cat/V:(abundances, times, volumes)}
            pickle.dump(results, file)
            file.flush()
    
except Exception as e:
    print(f"An error has occured: {e}")

finally: 
    print("Simulation run completed")
    
# --- BARRIDO DE K_0 ---

# Reaction constants:
k = [1]*8 # len(k)= # de reacciones
# Volume:
V = 1000
initial_food = [1000]*4 + [500]*2 + [0]*2 # initial molecules number
k_change = np.logspace(-5,-3,5)
print(k_change)
output_file = "barrido_k0_stream.pkl"

try:
    with open(output_file, "ab") as file:
        for k_i in k_change:
            k[0] = k_i
            abundances, times, volumes = chemistry(method, n_iterations, f,
                                    initial_food, k, V)
            print(f"Simulation completed for k_0 = {k_i}")
            results_k0 = {k_i : (abundances, times, volumes)}
            pickle.dump(results_k0, file)
            file.flush()
            
except Exception as e:
    print(f"An error has occured: {e}")

finally: 
    print("Simulation run completed")
# --- BARRIDO DE K_2 ---

# Reaction constants:
k = [1]*8 # len(k)= # de reacciones
# Volume:
V = 1000
initial_food = [1000]*4 + [500]*2 + [0]*2 # initial molecules number
final_volume = []
k_change = np.logspace(-5,-3,5)
print(k_change)
output_file = "barrido_k2_stream.pkl"

try: 
    with open(output_file, "ab") as file:

        for k_i in k_change:
            k[2] = k_i
            print(f"Starting simulation for k_2 = {k_i}")
            abundances, times, volumes = chemistry(method, n_iterations, f,
                                    initial_food, k, V)
            print(f"Simulation completed for k_2 = {k_i}")

            results_k2 = {k_i:(abundances, times, volumes)}
            pickle.dump(results_k2, file)

            file.flush()
except Exception as e:
    print(f"An error has occured: {e}")

finally: 
    print("Simulation run completed")


# -- BARRIDO RESTO DE K ---
k_change = np.logspace(-5,7,24)
for k_ns in[1,3,4,5,6,7]:
    print(f"STARTING SIMULATIONS FOR K_{k_ns}")
    barrido_k(k_ns, k_change)