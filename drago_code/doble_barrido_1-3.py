import numpy as np
import sys
import os
import sys

# Get the directory where your current script is located
current_dir = os.path.dirname(os.path.abspath(__file__))

# Move up one level to 'catalytic_world'
parent_dir = os.path.dirname(current_dir)

# Add it to the path
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from fun_gilles import *

import pickle

def barrido_k(k_n_1, k_n_2, k_values, num_replicates):
    # Reaction constants:
    k = [1]*8 # len(k)= # de reacciones
    # Volume:
    V = 100
    initial_food = [100]*4 + [50]*2 + [0]*2 # initial molecules number
    output_file = f"barrido_k{k_n_1}_k{k_n_2}.pkl"

    try: 
        with open(output_file, "ab") as file:
            i = 0
            for k_i in k_values:
                for k_j in k_values:
                    k[k_n_1] = k_i
                    k[k_n_2] = k_j
                    print(f"Starting simulation {i} for k_{k_n_1} = {k_i} and k_{k_n_2} = {k_j}")
                    replicate_results = []
                    
                    for r in range(num_replicates):
                        print(f"  -> Replicate {r+1}/{num_replicates}")
                        abundances, times, volumes = chemistry(method, n_iterations, f,
                                                                initial_food, k, V)
                        replicate_results.append((abundances, times, volumes))

                    # Store the list of replicates under the key k_i
                    results_ki = {i: replicate_results}
                    i += 1
                    print(f"Simulation {i} completed for k_{k_n_1} = {k_i} and k_{k_n_2} = {k_j}")

                    pickle.dump(results_ki, file)

                    file.flush()
    except Exception as e:
        print(f"An error has occured: {e}")

    finally: 
        print("Simulation run completed")
          

f = 'examples/reactions_autocat.txt'
reactions = read_file(f)
species = obtain_species(reactions)

n_reps = 1
n_iterations = 1e6
method = "Protocell" # Gillespie or Deterministic

# -- BARRIDO RESTO DE K ---
k = [1]*8 # len(k)= # de reacciones
# # Volume:
V = 100
initial_food = [100]*4 + [50]*2 + [0]*2 # initial molecules number
k_change = np.logspace(-5,7,13)
print(k_change)
k = [1]*8
barrido_k(1, 3, k_change, n_reps)