import numpy as np
import sys
sys.path.append('..')  # allow importing from parent directory
from fun_gilles import *
import pickle
import pandas as pd

def barrido_k(k_n, k_values, num_replicates):
    # Reaction constants:
    k = [1]*8 # len(k)= # de reacciones
    # Volume:
    V = 100
    initial_food = [100]*4 + [50]*2 + [0]*2 # initial molecules number
    output_file = f"barrido_k{k_n}_stream2.pkl"

    try: 
        with open(output_file, "ab") as file:

            for k_i in k_values:
                k[k_n] = k_i
                print(f"Starting simulation for k_{k_n} = {k_i}")
                replicate_results = []
                
                for r in range(num_replicates):
                    print(f"  -> Replicate {r+1}/{num_replicates}")
                    abundances, times, volumes = chemistry(method, n_iterations, f,
                                                            initial_food, k, V)
                    replicate_results.append((abundances, times, volumes))

                # Store the list of replicates under the key k_i
                results_ki = {k_i: replicate_results}
                
                print(f"Simulation completed for k_{k_n} = {k_i}")

                results_ki = {k_i:replicate_results}
                pickle.dump(results_ki, file)

                file.flush()
    except Exception as e:
        print(f"An error has occured: {e}")

    finally: 
        print("Simulation run completed")
        
def barrido_ab0(k_n, k_values, cat_abundance, num_replicates):
    # Reaction constants:
    k = [1]*8 # len(k)= # de reacciones
    # Volume:
    V = 100
    initial_food = [100]*4 + [0]*2 + [0]*2 # initial molecules number
    output_file = f"barrido_ab0_k{k_n}_{k_values}_stream2.pkl"
    k[k_n] = k_values
    try: 
        with open(output_file, "ab") as file:

            for cat in cat_abundance:
                initial_food[4] = cat
                print(f"Starting simulation for ab0 = {initial_food[4]/V}")
                replicate_results = []
                
                for r in range(num_replicates):
                    print(f"  -> Replicate {r+1}/{num_replicates}")
                    abundances, times, volumes = chemistry(method, n_iterations, f,
                                                            initial_food, k, V)
                    replicate_results.append((abundances, times, volumes))
                
                
                print(f"Simulation completed for ab0 = {initial_food[4]/V}")

                results = {initial_food[4]/V:replicate_results}
                pickle.dump(results, file)

                file.flush()
    except Exception as e:
        print(f"An error has occured: {e}")

    finally: 
        print("Simulation run completed")    

f = 'examples/reactions_autocat.txt'
reactions = read_file(f)
species = obtain_species(reactions)

n_reps = 10
# --- BARRIDO DE CONDICIONES INICIALES ---
n_iterations = 1e6
method = "Protocell" # Gillespie or Deterministic
# Reaction constants:
k = [1]*8 # len(k)= # de reacciones
# Volume:
V = 100
initial_food = [100]*4 + [0]*4 # initial molecules number

for k_1 in [1e-4,1,1e4]:
    cat_abundance= (np.logspace(-4, 6, 23) * 1000)[17:] 
    barrido_ab0(1, k_1, cat_abundance, n_reps)

for k_3 in [1e-4,1,1e4]:
    cat_abundance= np.logspace(-4, 6, 23) * 1000
    barrido_ab0(3, k_3, cat_abundance, n_reps)

# -- BARRIDO RESTO DE K ---
k = [1]*8 # len(k)= # de reacciones
# # Volume:
V = 100
initial_food = [100]*4 + [50]*2 + [0]*2 # initial molecules number
k_change = np.logspace(-5,7,25)
print(k_change)
for k_ns in range(0,8):
    k = [1]*8
    print(f"STARTING SIMULATIONS FOR K_{k_ns}")
    barrido_k(k_ns, k_change, n_reps)