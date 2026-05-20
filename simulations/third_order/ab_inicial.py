import os
import numpy as np
import sys
import pickle
from multiprocessing import Pool

sys.path.append('..')
from fun_gilles import *

# 1. Update helper to unpack 7 items
def run_single_simulation(args):
    # Now unpacking 7 items to match the tuple below
    method, n_iterations, f, initial_food, k, V, thresh = args
    return chemistry(method, n_iterations, f, initial_food, k, V, threshold=thresh)

def barrido_ab0(cat_abundance, num_replicates, processes):
    V = 100
    initial_food_base = [100]*4 + [0]*4
    output_file = "barrido_ab0_MM.pkl"
    k_vals = [1, 1, 1e2, 1]*2

    with Pool(processes=processes) as pool:
        with open(output_file, "ab") as file:
            for cat in cat_abundance:
                current_initial_food = initial_food_base.copy()
                current_initial_food[4] = cat
                
                # Correctly structured tuple (7 items)
                sim_args = (method, n_iterations, f, current_initial_food, k_vals, V, 0)
                task_list = [sim_args] * num_replicates

                results_list = []
                # imap_unordered is great for keeping memory low
                for result in pool.imap_unordered(run_single_simulation, task_list):
                    results_list.append(result)

                pickle.dump({cat/V: results_list}, file)
                file.flush()
                print(f"Done ab0 = {cat/V}")

# Global parameters
f = '../examples/reactions_autocat.txt'
n_reps = 1
n_iterations = 1e5
method = "Protocell"

if __name__ == "__main__":
    n_cpus = int(os.getenv('SLURM_CPUS_PER_TASK', 1))
    cat_abundance = np.linspace(0, 1000, 21) * 100
    barrido_ab0(cat_abundance, n_reps, n_cpus)