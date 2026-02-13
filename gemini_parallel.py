import numpy as np
import sys
import pickle
import os
from multiprocessing import Pool

# Ensure parent directory is accessible
sys.path.append('..')
from fun_gilles import *

def run_single_simulation(args):
    method, n_iterations, f, initial_food, k, V = args
    # Ensure threshold is handled
    return chemistry(method, n_iterations, f, initial_food, k, V)

def barrido_k(k_n, k_values, num_replicates, processes):
    V = 100
    initial_food = [100]*4 + [50]*2 + [0]*2
    output_file = f"barrido_k{k_n}_parallel.pkl"
    k_base = [1]*8

    # Open pool ONCE outside the loop to save memory and overhead
    with Pool(processes=processes) as pool:
        with open(output_file, "ab") as file:
            for k_i in k_values:
                k_base[k_n] = k_i
                print(f"Starting k_{k_n} = {k_i}")
                
                sim_args = (method, n_iterations, f, initial_food, k_base.copy(), V)
                task_list = [sim_args] * num_replicates

                # Use imap_unordered to save memory
                results_list = []
                for result in pool.imap_unordered(run_single_simulation, task_list):
                    results_list.append(result)
                
                pickle.dump({k_i: results_list}, file)
                file.flush()
                print(f"Completed k_{k_n} = {k_i}")

def barrido_ab0(k_n, k_val, cat_abundance, num_replicates, processes):
    V = 100
    initial_food_base = [100]*4 + [0]*2 + [0]*2
    output_file = f"barrido_ab0_k{k_n}_{k_val}_parallel.pkl"
    k = [1]*8
    k[k_n] = k_val

    with Pool(processes=processes) as pool:
        with open(output_file, "ab") as file:
            for cat in cat_abundance:
                current_initial_food = initial_food_base.copy()
                current_initial_food[4] = cat
                print(f"Starting ab0 = {cat/V}")
                
                sim_args = (method, n_iterations, f, current_initial_food, k, V)
                task_list = [sim_args] * num_replicates

                results_list = []
                for result in pool.imap_unordered(run_single_simulation, task_list):
                    results_list.append(result)

                pickle.dump({cat/V: results_list}, file)
                file.flush()

# Global parameters
f = '../examples/reactions_autocat.txt'
n_reps = 10
n_iterations = 1e7
method = "Protocell"

if __name__ == "__main__":
    # Get CPUs from Slurm, default to 1 if running locally
    n_cpus = int(os.getenv('SLURM_CPUS_PER_TASK', 1))
    
    # --- BARRIDO DE CONDICIONES INICIALES ---
    for k_1 in [1e-4]:
        cat_abundance = np.linspace(0, 100, 21)* 100
        barrido_ab0(1, k_1, cat_abundance, n_reps, n_cpus)

    # -- BARRIDO RESTO DE K ---
    k_change = np.logspace(-5, 7, 25)
    for k_ns in range(0, 8):
        print(f"STARTING SIMULATIONS FOR K_{k_ns}")
        barrido_k(k_ns, k_change, n_reps, n_cpus)