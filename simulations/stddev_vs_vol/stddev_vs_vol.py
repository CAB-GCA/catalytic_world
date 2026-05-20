import numpy as np
import sys
import pickle
import os
sys.path.append('..')
from fun_gilles import *

# --- 1. CONFIGURATION ---
VOLUME_EXPONENTS = np.arange(2, 5.5, 0.5)
VOLUMES = 10**VOLUME_EXPONENTS
REACTIONS_FILE = "../examples/reactions_autocat.txt" # Ensure this is in the same folder or use absolute path
K_CONSTANTS = [1, 1, 1, 1, 1, 1, 1, 1] 
INITIAL_CONCENTRATIONS = {
    'a': 1.0, 'b': 1.0, 'c': 1.0, 'd': 1.0,
    'ab': 0.30598098, 'cd': 0.30856037, 
    'cab': 0.19095599, 'acd': 0.19450266
}

# --- 2. SLURM ARRAY LOGIC ---
try:
    # SLURM_ARRAY_TASK_ID is the standard way to get the index on Drago
    task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", sys.argv[1] if len(sys.argv) > 1 else 0))
    V = VOLUMES[task_id]
except (IndexError, ValueError):
    print("Error: Invalid Task ID")
    sys.exit(1)

V_int = int(round(V))
MAX_ITERATIONS = int(V_int * 20) # Increased to ensure steady state is reached
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
        # We set threshold=0 to force the simulation to run for the full iterations
        # or until your internal counter triggers.
        abundances, times, final_V = chemistry(
            method='Protocell', 
            iterations=MAX_ITERATIONS,
            file=REACTIONS_FILE, 
            initial_food=initial_abundances, 
            k=K_CONSTANTS, 
            V=V_int, 
            threshold=0 
        )
        
        if len(abundances) >= BLOCK_SIZE:
            # Concentration calculation: [Conc] = Molecules / Volume
            # final_V is a 1D array of volume over time
            abundances_slice = abundances[-BLOCK_SIZE:, :]
            V_slice = final_V[-BLOCK_SIZE:]
            
            # Divide each row of abundances by the corresponding volume at that time
            last_block_concentrations = (abundances_slice.T / V_slice).T
            
            # Fixed call to match fun_gilles.py signature: block_statistics(abundances)
            block_std = block_statistics(last_block_concentrations)
            results_for_this_volume.append(block_std.tolist())
            print(f"Rep {r} completed.")
            
    except Exception as e:
        print(f"Rep {r} failed: {e}")

# --- 4. SAVE ---
output_filename = f"resultsv5_V_{V_int}.pkl"
with open(output_filename, "wb") as f:
    pickle.dump({V_int: results_for_this_volume}, f)

print(f"Finished. Data saved to {output_filename}")