###############################################################################################
# FOR AB0_SWEEP THE WORKFLOW IS THE FOLLOWING:
    # 1. GENERATE THE PKL FILES USING ab0_sweep_pkl_generator.py SETTING THE DESIRED PARAMETERS
    # 2. RUN ab0_sweep.py TO GET THE FIGURES AND STATS
    # 3. TO CHANGE THE TITLES OF THE FIGURE (WHICH MAY NOT CORRESPOND TO THE CHOSEN 
    #                                        PARAMETERS) GO TO ../src/plotting.py
###############################################################################################
import numpy as np
import pickle
import multiprocessing as mp
import os
import sys
import time

# --- ROBUST IMPORT FIX ---
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from src.fun_gilles import *

def run_single_simulation(args):
    """
    Wrapper function to execute a single Gillespie run.
    'args' is a tuple containing all the necessary parameters for the simulation.
    """
    # Unpack arguments
    (condition_val, rep_id, abundances_init, m, k_types, k, c, 
     V_init, iterations, c_reactants, threshold) = args
    
    # We update the specific condition we are sweeping. 
    # Assuming AB is at index 0 (UPDATE THIS INDEX if AB is elsewhere in your array!)
    abundances_init_copy = abundances_init.copy()
    abundances_init_copy[0, 4] = condition_val * V_init # Assuming condition is concentration
    
    # Run the optimized Gillespie function
    try:
        a, t, v = gillespieProtocell(
            abundances_init_copy, m, k_types, k, c, V_init, iterations, c_reactants, threshold
        )
        return (condition_val, rep_id, (a, t, v))
    except Exception as e:
        print(f"Simulation failed for condition {condition_val}, rep {rep_id}: {e}")
        return (condition_val, rep_id, None)

def main():
    # --- 1. SETUP SIMULATION PARAMETERS ---
    file_path = os.path.join(project_root, 'examples', 'reactions_autocat.txt')
    reactions = read_file(file_path)
    header = get_header(file_path)
    species = obtain_species(reactions)
    c = c_matrix(reactions, species, header)
    c_reactants = reactants(reactions, species, header)
    k_types = reactions[:, -1]
    m = reactions.shape[0]
    
    # Base kinetic constants (Update these to match the specific k1 sweep you are doing)
    # e.g., k1 = 1.0
    k = [1,1,1,1]+[1,1,1,1] 
    
    V_init = 1000.0
    iterations = 500000
    
    # Initialize food and standard species
    abundances_init = np.zeros((1, len(species)))
    abundances_init[0, :4] = 1 * V_init 
    
    # --- 2. DEFINE THE SWEEP ---
    # Example conditions matching your plots (0.0 to 10.0)
    ab0_conditions = np.linspace(0.0, 10.0, 21) 
    replicates_per_condition = 10
    output_filename = "barrido_ab0_k1_1_local.pkl"
    threshold = threshold_function(V_init, 0.05) 
    # Prepare all tasks for the CPU pool
    tasks = []
    for cond in ab0_conditions:
        for rep in range(replicates_per_condition):
            task = (cond, rep, abundances_init, m, k_types, k, c, 
                    V_init, iterations, c_reactants, threshold)
            tasks.append(task)
            
    print(f"Total simulations to run: {len(tasks)}")
    
    # --- 3. EXECUTE IN PARALLEL ---
    # Use max CPUs minus 1 so your computer doesn't completely freeze
    num_cores = max(1, mp.cpu_count() - 1) 
    print(f"Firing up {num_cores} CPU cores...\n")
    
    results_dict = {cond: [] for cond in ab0_conditions}
    
    start_time = time.time()
    
    # Use standard Pool to process tasks
    with mp.Pool(processes=num_cores) as pool:
        for i, result in enumerate(pool.imap_unordered(run_single_simulation, tasks)):
            cond, rep_id, data_tuple = result
            
            if data_tuple is not None:
                results_dict[cond].append(data_tuple)
                
            # Simple progress tracker
            if (i + 1) % 10 == 0:
                print(f"Completed {i + 1}/{len(tasks)} runs...")

    # --- 4. SAVE RESULTS ---
    output_path = os.path.join(project_root, 'data_archive', output_filename)
    with open(output_path, 'wb') as f:
        pickle.dump(results_dict, f)
        
    elapsed = (time.time() - start_time) / 60
    print(f"\nAll done! Saved to {output_path} in {elapsed:.2f} minutes.")

if __name__ == '__main__':
    # Required for Windows multiprocessing
    mp.freeze_support()
    main()