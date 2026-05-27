#  FOR THE VISUALIZATION RUN "k_sweep_2d.py"

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
        
    start_idx = round(len(volumes) / 2) if values_used is None else 0
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
    # iterations = 1000
    iterations = int(2e5)
    threshold = 0
    
    # Standard initial conditions for autocat
    abundances_init = np.zeros((1, len(species)))
    abundances_init[0, :4] = 1.0 * V_init 
    abundances_init[0, 4] = 1.0 * V_init # ab0 = 1
    
    # --- DEFINE THE 2D SWEEP ---
    k_change = np.logspace(-5, 7, 25)
    replicates_per_condition = 2
    
    # Target indices in the 'k' array
    # To sweep k2 vs k4 (which are indices 0 and 2 in Python)
    # k_y_index = 0
    # k_x_index = 2 
    # To sweep k2 vs k4 (which are indices 1 and 3 in Python)
    k_y_index = 1 
    k_x_index = 3 
    
    prefix = f"sweep_k{k_y_index}_k{k_x_index}"
    
    # Prepare task list
    tasks = []
    for row_idx, k_y_val in enumerate(k_change):
        for col_idx, k_x_val in enumerate(k_change):
            for rep in range(replicates_per_condition):
                task = (row_idx, col_idx, rep, k_y_val, k_x_val, k_y_index, k_x_index, 
                        abundances_init, m, k_types, base_k, c, V_init, iterations, c_reactants, threshold)
                tasks.append(task)
                
    total_tasks = len(tasks)
    print(f"Total simulations to run: {total_tasks}")
    
    # Initialize grids to store the sum of alpha/R2 (we will average them later)
    grid_shape = (len(k_change), len(k_change))
    alpha_sum_grid = np.zeros(grid_shape)
    r2_sum_grid = np.zeros(grid_shape)
    counts_grid = np.zeros(grid_shape)
    
    # --- EXECUTE IN PARALLEL ---
    num_cores = max(1, mp.cpu_count() - 1)
    print(f"Firing up {num_cores} CPU cores...")
    start_time = time.time()
    
    with mp.Pool(processes=num_cores) as pool:
        for i, result in enumerate(pool.imap_unordered(run_2d_simulation, tasks)):
            row_idx, col_idx, rep_id, alpha, r2 = result
            
            # Accumulate results
            alpha_sum_grid[row_idx, col_idx] += alpha
            r2_sum_grid[row_idx, col_idx] += r2
            counts_grid[row_idx, col_idx] += 1
            
            if (i + 1) % 100 == 0:
                percent = ((i+1) / total_tasks) * 100
                print(f"Completed {i + 1}/{total_tasks} runs ({percent:.1f}%)...")

    # --- CALCULATE MEANS & SAVE ---
    print("\nCalculating final grids...")
    # Avoid division by zero just in case
    counts_grid[counts_grid == 0] = 1 
    
    alpha_mean_grid = alpha_sum_grid / counts_grid
    r2_mean_grid = r2_sum_grid / counts_grid
    
    out_dir = os.path.join(project_root, 'data_archive')
    os.makedirs(out_dir, exist_ok=True)
    
    alpha_file = os.path.join(out_dir, f"alpha_{prefix}.pkl")
    r2_file = os.path.join(out_dir, f"r2_{prefix}.pkl")
    
    with open(alpha_file, "wb") as f:
        pickle.dump(alpha_mean_grid, f)
        
    with open(r2_file, "wb") as f:
        pickle.dump(r2_mean_grid, f)

    elapsed = (time.time() - start_time) / 60
    print(f"\nAll done in {elapsed:.2f} minutes!")
    print(f"Saved to:\n - {alpha_file}\n - {r2_file}")


if __name__ == '__main__':
    mp.freeze_support()
    main()