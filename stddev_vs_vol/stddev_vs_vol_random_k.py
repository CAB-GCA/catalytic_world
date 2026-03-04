import numpy as np
import sys
import pickle
import os

sys.path.append('..')
from fun_gilles import *

# --- 1. CONFIGURATION ---
VOLUME_EXPONENTS = np.arange(2, 5.5, 0.5) # 7 volumes
VOLUMES = 10**VOLUME_EXPONENTS
REACTIONS_FILE = "../examples/reactions_autocat.txt"
INITIAL_CONCENTRATIONS = {
    'a': 1.0, 'b': 1.0, 'c': 1.0, 'd': 1.0,
    'ab': 0.30598098, 'cd': 0.30856037, 
    'cab': 0.19095599, 'acd': 0.19450266
}

N_K_SETS = 10  # How many different random sets of k constants to try
REPS = 15      # Replicates per (Volume, K-set)
BLOCK_SIZE = 500

# --- 2. HPC ARRAY LOGIC ---
# Total tasks should be len(VOLUMES) * N_K_SETS (e.g., 7 * 10 = 70 tasks)
try:
    task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", sys.argv[1] if len(sys.argv) > 1 else 0))
    
    # Map task_id to Volume index and K-set index
    v_idx = task_id % len(VOLUMES)
    k_idx = (task_id // len(VOLUMES))    
    V = VOLUMES[v_idx]
except (IndexError, ValueError):
    print(f"Error: Invalid Task ID {task_id}")
    sys.exit(1)

V_int = int(round(V))
MAX_ITERATIONS = int(V_int * 20) 

# --- 3. RANDOM K GENERATION ---
# Seed the RNG with k_idx so all volumes in a "set" use the same K but different sets are unique
np.random.seed(k_idx)

# Initialize an array of 8 zeros
K_CONSTANTS = np.zeros(8)

# Indices: 0, 2, 4, 6 (Evens) -> Higher than 1e0 (10^0 to 10^7)
even_exponents = np.random.uniform(0, 7, size=4)
K_CONSTANTS[0::2] = 10**even_exponents

# Indices: 1, 3, 5, 7 (Odds) -> Lower than 1e0 (10^-7 to 10^0)
odd_exponents = np.random.uniform(-7, 0, size=4)
K_CONSTANTS[1::2] = 10**odd_exponents

# Get species ordering
reactions = read_file(REACTIONS_FILE)
SPECIES = obtain_species(reactions)

# --- 4. SIMULATION ---
results_for_this_run = []
initial_abundances = np.array([round(INITIAL_CONCENTRATIONS.get(s, 0.0) * V_int) for s in SPECIES])

print(f"Task {task_id}: Volume {V_int}, K-Set {k_idx}")
print(f"Rates used: {K_CONSTANTS.round(3)}")

for r in range(REPS):
    try:
        abundances, times, final_V = chemistry(
            method='Protocell', 
            iterations=max(MAX_ITERATIONS, BLOCK_SIZE + 100),
            file=REACTIONS_FILE, 
            initial_food=initial_abundances, 
            k=K_CONSTANTS.tolist(), 
            V=V_int, 
            threshold=0
        )
        
        if len(abundances) >= BLOCK_SIZE:
            abundances_slice = abundances[-BLOCK_SIZE:, :]
            V_slice = final_V[-BLOCK_SIZE:]
            last_block_concentrations = abundances_slice / V_slice[:, np.newaxis]
            
            # Using standard signature from fun_gilles.py
            block_std = block_statistics(last_block_concentrations)
            results_for_this_run.append(block_std.tolist())
            
    except Exception as e:
        print(f"Rep {r} failed: {e}")

# --- 5. SAVE ---
# Save with both V and K-index to avoid collisions
output_filename = f"results_random_V{V_int}_K{k_idx}.pkl"
with open(output_filename, "wb") as f:
    pickle.dump({
        'volume': V_int,
        'k_idx': k_idx,
        'k_values': K_CONSTANTS,
        'stds': results_for_this_run
    }, f)

print(f"Finished. Data saved to {output_filename}")
