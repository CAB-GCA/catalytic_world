import numpy as np
import sys
import pickle
import os

sys.path.append('..')
from fun_gilles import *

# --- 1. CONFIGURATION ---
VOLUME_EXPONENTS = np.arange(2, 5.5, 0.5)
VOLUMES = 10**VOLUME_EXPONENTS
REACTIONS_FILE = "../examples/reactions_autocat.txt"
K_CONSTANTS = [1.0] * 8
INITIAL_CONCENTRATIONS = {
    'a': 1.0, 
    'b': 1.0, 
    'c': 1.0, 
    'd': 1.0,
    'ab': 0.30598098, 
    'cd': 0.30856037, 
    'cab': 0.19095599, 
    'acd': 0.19450266
}

# --- 2. HPC ARRAY LOGIC ---
try:
    task_id = int(sys.argv[1])  # Slurm passes this
    V = VOLUMES[task_id]
except (IndexError, ValueError):
    print("Error: Provide a valid index (0 to 6) as an argument.")
    sys.exit(1)

V_int = int(round(V))
MAX_ITERATIONS = int(V_int * 10)
REPS = 5
BLOCK_SIZE = 500

# Get species ordering
reactions = read_file(REACTIONS_FILE)
SPECIES = obtain_species(reactions)

# --- 3. SIMULATION ---
results_for_this_volume = []
initial_abundances = np.array([round(INITIAL_CONCENTRATIONS.get(s, 0.0) * V_int) for s in SPECIES])

print(f"Starting Volume {V_int} (Task {task_id})")

for r in range(REPS):
    try:
        abundances, times, final_V = chemistry(
            method='Protocell', 
            iterations=max(MAX_ITERATIONS, BLOCK_SIZE + 100),
            file=REACTIONS_FILE, 
            initial_food=initial_abundances, 
            k=K_CONSTANTS, 
            V=V_int, 
            threshold=0
        )
        
        if len(abundances) >= BLOCK_SIZE:
            # Concentration calculation
            abundances_slice = abundances[-BLOCK_SIZE:, :]
            V_slice = final_V[-BLOCK_SIZE:]
            last_block_concentrations = abundances_slice / V_slice[:, np.newaxis]
            
            # Use the reshape-based block_statistics we discussed
            block_std = block_statistics(last_block_concentrations, axis= 0)
            results_for_this_volume.append(block_std.tolist())
            
    except Exception as e:
        print(f"Rep {r} failed: {e}")

# --- 4. SAVE UNIQUE FILE ---
# Save as result_V_100.pkl etc. to avoid file writing conflicts
output_filename = f"results_task_{task_id}_V_{V_int}.pkl"

with open(output_filename, "wb") as f:
    pickle.dump({V_int: results_for_this_volume}, f)

print(f"Finished. Data saved to {output_filename}")