import pickle
import sys
import os
import numpy as np

# --- Path setup ---
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from fun_gilles import *

def xyeq_vs_c_ar(initial_c, initial_food, method, n_iterations, 
                 k, V0, k_var, n_reps=5):
    
    f_path = "../examples/reactions_XYC.txt"

    for j in range(len(k_var)):
        current_kar = k_var[j]
        k[1] = current_kar 
        
        # Matrix to store ALL results: rows = concentrations, columns = replicates
        # Shape: (len(initial_c), n_reps)
        xy_results_matrix = np.zeros((len(initial_c), n_reps))
        
        for i in range(len(initial_c)):
            initial_food[0] = initial_c[i]
            
            for r in range(n_reps):
                abundances, times, V = chemistry(method, n_iterations, f_path, 
                                                initial_food, k, V0)
                # Store the raw result for this specific replicate
                xy_results_matrix[i, r] = abundances[-1, -1] / V0
            
            if i % 10 == 0:
                print(f"k_ar: {current_kar:.1e} | Progress: {i}/{len(initial_c)}")

        # Save the full matrix for this k_ar
        with open("XYC_gillespie_kar.pkl", "ab") as file:
            pickle.dump({current_kar : xy_results_matrix}, file)
            file.flush()
        
        print(f"--- Finished all replicates for k_a_r = {current_kar} ---")

# --- Parameters ---
k = [1, 1, 1, 1] 
V0 = 1000
initial_food = [0, 1000, 1000] + [0]*2 
initial_c = np.linspace(0, 5000, 51)
initial_c[0] = 1                     # Avoid c=0 because then all probabilities are zero and the simulation can't proceed
k_var = np.logspace(-4, 3, 8)        
n_iterations = int(1e5)
n_reps = 10                          # Number of replicates per point
method = "Gillespie"

# Execute
xyeq_vs_c_ar(initial_c, initial_food, method, n_iterations,
             k, V0, k_var, n_reps=n_reps)

print(f"\nSimulation Complete.\nReplicates per point: {n_reps}")