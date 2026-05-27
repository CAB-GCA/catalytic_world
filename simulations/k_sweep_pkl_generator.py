# TO VISUALIZE THE RESULTS RUN "k_sweep.py"
import numpy as np
import pickle
import multiprocessing as mp
import os
import sys
import time

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from src.fun_gilles import *

def run_single_simulation(args):
    """Worker: Runs a single Gillespie simulation for a specific k value."""
    (k_val, rep_id, target_k_index, abundances_init, m, k_types, base_k, c, 
     V_init, iterations, c_reactants, threshold) = args
    
    # Inject the specific sweep condition into the kinetic constants array
    k_local = base_k.copy()
    k_local[target_k_index] = k_val
    
    try:
        a, t, v = gillespieProtocell(
            abundances_init, m, k_types, k_local, c, V_init, iterations, c_reactants, threshold
        )
        # Return the condition, replicate ID, and the data tuple
        return (k_val, rep_id, (a, t, v))
    except Exception as e:
        print(f"Simulation failed for k_{target_k_index+1} = {k_val}, rep {rep_id}: {e}")
        return (k_val, rep_id, None)

def run_sweep_for_index(target_k_index, k_conditions, replicates_per_condition, base_params):
    """Manages the multiprocessing pool for a specific kinetic parameter."""
    abundances_init, m, k_types, base_k, c, V_init, iterations, c_reactants, threshold = base_params
    
    filename = f"barrido_k{target_k_index}.pkl"
    output_path = os.path.join(project_root, 'data_archive', filename)
    
    print(f"\n{'='*50}")
    print(f"Starting sweep for k_{target_k_index + 1}...")
    print(f"Target Output: {filename}")
    print(f"{'='*50}")
    
    tasks = []
    for k_val in k_conditions:
        for rep in range(replicates_per_condition):
            task = (k_val, rep, target_k_index, abundances_init, m, k_types, base_k, c, 
                    V_init, iterations, c_reactants, threshold)
            tasks.append(task)
            
    total_tasks = len(tasks)
    num_cores = max(1, mp.cpu_count() - 1)
    print(f"Total simulations: {total_tasks}. Firing up {num_cores} cores...")
    
    results_dict = {cond: [] for cond in k_conditions}
    start_time = time.time()
    
    with mp.Pool(processes=num_cores) as pool:
        for i, result in enumerate(pool.imap_unordered(run_single_simulation, tasks)):
            k_val, rep_id, data_tuple = result
            
            if data_tuple is not None:
                results_dict[k_val].append(data_tuple)
                
            if (i + 1) % 50 == 0:
                print(f"Completed {i + 1}/{total_tasks} runs...")

    # Save to disk using a streaming approach (chunk by condition) to mimic the cluster
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'wb') as f:
        for condition, replicates in results_dict.items():
            # Dump one condition dictionary at a time
            pickle.dump({condition: replicates}, f)
            
    elapsed = (time.time() - start_time) / 60
    print(f"Sweep for k_{target_k_index + 1} finished in {elapsed:.2f} minutes.")

def main():
    # --- 1. SETUP BASE SIMULATION PARAMETERS ---
    file_path = os.path.join(project_root, 'examples', 'reactions_autocat.txt')
    reactions = read_file(file_path)
    header = get_header(file_path)
    species = obtain_species(reactions)
    c = c_matrix(reactions, species, header)
    c_reactants = reactants(reactions, species, header)
    k_types = reactions[:, -1]
    m = reactions.shape[0]
    
    base_k = [1.0] * 8
    V_init = 1000.0
    iterations = 100000
    threshold = 0
    
    # Initial conditions for standard autocat
    abundances_init = np.zeros((1, len(species)))
    abundances_init[0, :4] = 1.0 * V_init  # A, B, C, D
    abundances_init[0, 4] = 1.0 * V_init   # AB
    
    base_params = (abundances_init, m, k_types, base_k, c, V_init, iterations, c_reactants, threshold)
    
    # --- 2. DEFINE THE SWEEPS ---
    replicates_per_condition = 10
    
    # You can adjust these logspaces based on your exact desired limits!
    # Sweep for k1 (Index 0): From 10^-4 to 10^7
    k0_conditions = np.logspace(-4, 7, 23)
    
    # Sweep for k2 (Index 1): From 10^-4 to 10^7
    k1_conditions = np.logspace(-4, 7, 23)
    
    # --- 3. EXECUTE ---
    run_sweep_for_index(target_k_index=0, k_conditions=k0_conditions, 
                        replicates_per_condition=replicates_per_condition, base_params=base_params)
                        
    run_sweep_for_index(target_k_index=1, k_conditions=k1_conditions, 
                        replicates_per_condition=replicates_per_condition, base_params=base_params)
                        
    print("\nAll 1D parameter sweeps generated successfully!")

if __name__ == '__main__':
    mp.freeze_support()
    main()