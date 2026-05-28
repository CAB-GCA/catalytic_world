# TO GENERATE THE PKL FILES RUN "k_sweep_2d_pkl_generator.py"

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

# --- 1. WORKER FUNCTIONS ---

def _get_alpha_from_vectors(times, volumes, values_used= None):
    """Calculates alpha and R^2 directly from vectors to save memory."""
    if len(volumes) < 10:
        return 0.0, 0.0
        
    start_idx = round(len(volumes) / 5) if values_used is None else 0
    end_idx = round(len(volumes) * values_used) if values_used is not None else len(volumes)
    
    log_v = np.log(volumes[start_idx:end_idx])
    t_subset = times[start_idx:end_idx]
    
    if len(t_subset) < 2:
        return 0.0, 0.0
        
    coefficients = np.polyfit(t_subset, log_v, 1)
    m_fit = coefficients[0]
    log_A_fit = coefficients[1]
    
    log_v_predicted = log_A_fit + m_fit * t_subset
    SS_res = np.sum((log_v - log_v_predicted)**2)
    SS_tot = np.sum((log_v - np.mean(log_v))**2)
    R_squared = 1 - (SS_res / SS_tot) if SS_tot != 0 else 0
    
    return m_fit, R_squared

def run_2d_simulation(args):
    """Worker: Runs simulation and IMMEDIATELY extracts Alpha/R2 to save RAM."""
    (row_idx, col_idx, rep_id, k_y_val, k_x_val, k_y_index, k_x_index, 
     abundances_init, m, k_types, base_k, c, V_init, iterations, c_reactants, threshold) = args
    
    # Create a fresh copy of the kinetic constants and inject the sweep values
    k_local = base_k.copy()
    k_local[k_y_index] = k_y_val
    k_local[k_x_index] = k_x_val
    
    try:
        a, t, v = gillespieProtocell(
            abundances_init, m, k_types, k_local, c, V_init, iterations, c_reactants, threshold
        )
        alpha, r2 = _get_alpha_from_vectors(t, v, values_used=1)
        
        return (row_idx, col_idx, rep_id, alpha, r2)
        
    except Exception as e:
        # If it fails (e.g. probability 0), return 0 alpha and 0 R2
        return (row_idx, col_idx, rep_id, 0.0, 0.0)


# --- 2. MAIN EXECUTION ---

def main():
    print("Initializing 2D Parameter Sweep...")
    
    # --- SETUP SIMULATION PARAMETERS ---
    file_path = os.path.join(project_root, 'examples', 'reactions_autocat.txt')
    reactions = read_file(file_path)
    header = get_header(file_path)
    species = obtain_species(reactions)
    c = c_matrix(reactions, species, header)
    c_reactants = reactants(reactions, species, header)
    k_types = reactions[:, -1]
    m = reactions.shape[0]
    
    # Base parameters
    base_k = [1.0] * 8
    V_init = 1000.0
    iterations = int(1e5)
    threshold = None
    
    # Standard initial conditions for autocat
    abundances_init = np.zeros((1, len(species)))
    abundances_init[0, :4] = 1.0 * V_init 
    abundances_init[0, 4] = 1.0 * V_init # ab0 = 1
    
    # --- DEFINE THE 2D SWEEP ---
    k_change = np.logspace(-5, 7, 25)
    replicates_per_condition = 5
    
    # Target indices in the 'k' array 
    # for k2 and k4 keep the following lines
    k_y_index = 1 
    k_x_index = 3 
    # for k1 and k3 keep the following lines
    # k_y_index = 0
    # k_x_index = 2    
    prefix = f"sweep_k{k_y_index}_k{k_x_index}"
    out_dir = os.path.join(project_root, 'data_archive')
    os.makedirs(out_dir, exist_ok=True)
    
    checkpoint_file = os.path.join(out_dir, f"checkpoint_{prefix}.pkl")
    
    # --- LOAD CHECKPOINT IF IT EXISTS ---
    if os.path.exists(checkpoint_file):
        print(f"Found existing checkpoint file! Loading previously completed simulations...")
        with open(checkpoint_file, "rb") as f:
            completed_tasks = pickle.load(f)
    else:
        completed_tasks = {}

    # --- PREPARE TASK LIST (Skipping completed ones) ---
    tasks = []
    for row_idx, k_y_val in enumerate(k_change):
        for col_idx, k_x_val in enumerate(k_change):
            for rep in range(replicates_per_condition):
                task_id = (row_idx, col_idx, rep)
                
                # Only add to tasks if it hasn't been done yet
                if task_id not in completed_tasks:
                    task = (row_idx, col_idx, rep, k_y_val, k_x_val, k_y_index, k_x_index, 
                            abundances_init, m, k_types, base_k, c, V_init, iterations, c_reactants, threshold)
                    tasks.append(task)
                
    total_new_tasks = len(tasks)
    print(f"Total simulations remaining to run: {total_new_tasks}")
    
    if total_new_tasks > 0:
        # --- EXECUTE IN PARALLEL (Optimized for OS responsiveness) ---
        # Leave at least 2 cores free for your operating system so it doesn't freeze.
        num_cores = max(1, mp.cpu_count() - 2) 
        print(f"Firing up {num_cores} CPU cores ")
        start_time = time.time()
        
        with mp.Pool(processes=num_cores) as pool:
            # chunksize=10 tells Python to hand out jobs in batches rather than 1 by 1.
            # This drastically reduces overhead and stops your computer from choking on communication.
            for i, result in enumerate(pool.imap_unordered(run_2d_simulation, tasks, chunksize=10)):
                row_idx, col_idx, rep_id, alpha, r2 = result
                
                # Save to our dictionary
                completed_tasks[(row_idx, col_idx, rep_id)] = (alpha, r2)
                
                # Save the checkpoint to disk every 100 runs so we don't lose progress if it crashes
                if (i + 1) % 100 == 0:
                    percent = ((i+1) / total_new_tasks) * 100
                    print(f"Completed {i + 1}/{total_new_tasks} remaining runs ({percent:.1f}%)...")
                    with open(checkpoint_file, "wb") as f:
                        pickle.dump(completed_tasks, f)
        
        # Final checkpoint save when all tasks are done
        with open(checkpoint_file, "wb") as f:
            pickle.dump(completed_tasks, f)
            
        elapsed = (time.time() - start_time) / 60
        print(f"\nFinished running new simulations in {elapsed:.2f} minutes!")
    else:
        print("All simulations were already completed! Proceeding to generate final grids...")

    # --- RECONSTRUCT FINAL GRIDS FROM THE CHECKPOINT ---
    print("\nCalculating final grids...")
    grid_shape = (len(k_change), len(k_change))
    alpha_sum_grid = np.zeros(grid_shape)
    r2_sum_grid = np.zeros(grid_shape)
    counts_grid = np.zeros(grid_shape)

    for (row, col, rep), (alpha, r2) in completed_tasks.items():
        alpha_sum_grid[row, col] += alpha
        r2_sum_grid[row, col] += r2
        counts_grid[row, col] += 1

    # Avoid division by zero
    counts_grid[counts_grid == 0] = 1 
    
    alpha_mean_grid = alpha_sum_grid / counts_grid
    r2_mean_grid = r2_sum_grid / counts_grid
    
    alpha_file = os.path.join(out_dir, f"alpha_{prefix}.pkl")
    r2_file = os.path.join(out_dir, f"r2_{prefix}.pkl")
    
    with open(alpha_file, "wb") as f:
        pickle.dump(alpha_mean_grid, f)
        
    with open(r2_file, "wb") as f:
        pickle.dump(r2_mean_grid, f)

    print(f"Saved to:\n - {alpha_file}\n - {r2_file}")

if __name__ == '__main__':
    mp.freeze_support()
    main()