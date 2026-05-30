# TEST CODE TO SEE RESULTS USING THE DETERMINISTIC APPROACH
# TO GENERATE THE PKL FILES RUN "k_sweep_2d_pkl_generator.py"

import numpy as np
import pickle
import multiprocessing as mp
import os
import time
import warnings
from scipy.integrate import solve_ivp

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)

# --- 1. ODE SYSTEM ---

def odes_intermediates(t, X, K_intermediates, fi):
    """Deterministic system of ODEs for the minimal autocatalytic cycle."""
    V, ab, cd, acd, cab = X
    
    # Food concentrations (constant)
    a, b, c, d = 1.0, 1.0, 1.0, 1.0
    k1, k2, k3, k4, k5, k6, k7, k8 = K_intermediates
    
    # Volume derivative
    dVdt = (V / fi) * (k3 * b * acd + k7 * d * cab - k4 * ab * cd - k8 * ab * cd)
    
    # Species derivatives
    dabdt = k3 * b * acd + k6 * cab + k7 * cab * d - k4 * ab * cd - k5 * ab * c - k8 * cd * ab - (ab / V) * dVdt
    dcddt = k3 * b * acd + k7 * cab * d + k2 * acd - k4 * ab * cd - k1 * a * cd - k8 * cd * ab - (cd / V) * dVdt
    dacddt = k4 * ab * cd + k1 * a * cd - k3 * b * acd - k2 * acd - (acd / V) * dVdt
    dcabdt = k5 * ab * c + k8 * cd * ab - k6 * cab - k7 * cab * d - (cab / V) * dVdt
    
    return [dVdt, dabdt, dcddt, dacddt, dcabdt]


# --- 2. WORKER FUNCTIONS ---

def _get_alpha_from_vectors(times, volumes, values_used=None):
    """Calculates alpha and R^2 directly from vectors to save memory."""
    if len(volumes) < 10 or np.any(volumes <= 0):
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
    """Worker: Runs deterministic simulation and extracts Alpha/R2."""
    (row_idx, col_idx, k_y_val, k_x_val, k_y_index, k_x_index, 
     X_init, base_k, fi, t_span, t_eval) = args
    
    # Create a fresh copy of the kinetic constants and inject the sweep values
    k_local = base_k.copy()
    k_local[k_y_index] = k_y_val
    k_local[k_x_index] = k_x_val
    
    try:
        # Ignore RuntimeWarnings for zero division or overflow during ODE solving
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            # Solve the ODE using the BDF method (best for stiff chemical systems)
            res = solve_ivp(
                odes_intermediates, 
                t_span, 
                X_init, 
                method='BDF', 
                t_eval=t_eval, 
                args=(k_local, fi)
            )
            
        if not res.success:
            return (row_idx, col_idx, 0.0, 0.0)

        t = res.t
        v = res.y[0] # Volume is the first variable in the X array
        
        alpha, r2 = _get_alpha_from_vectors(t, v, values_used=1)
        
        return (row_idx, col_idx, alpha, r2)
        
    except Exception as e:
        # If integration fails, return 0
        return (row_idx, col_idx, 0.0, 0.0)


# --- 3. MAIN EXECUTION ---

def main():
    print("Initializing Deterministic 2D Parameter Sweep...")
    
    # --- SETUP SIMULATION PARAMETERS ---
    # Initial conditions from TFM (V0=100, [ab]0=0.5, [cd]0=0.5)
    V_init = 100.0
    ab_init = 0.5
    cd_init = 0.5
    acd_init = 0.0
    cab_init = 0.0
    
    X_init = [V_init, ab_init, cd_init, acd_init, cab_init]
    
    # fi calculation: Sum of initial concentrations
    fi = ab_init + cd_init + acd_init + cab_init
    
    # Base parameters (k1 through k8)
    base_k = [1.0] * 8
    
    # Time limits for integration
    t_span = (0, 500)
    t_eval = np.linspace(t_span[0], t_span[1], 1000)
    
    # --- DEFINE THE 2D SWEEP ---
    k_change = np.logspace(-5, 7, 25)
    
    # Target indices in the 'k' array 
    # For k2 and k4 (index 1 and 3)
    k_y_index = 1 
    k_x_index = 3 
    
    # For k1 and k3 uncomment below:
    # k_y_index = 0
    # k_x_index = 2    
    
    prefix = f"sweep_k{k_y_index}_k{k_x_index}_deterministic"
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
            
            task_id = (row_idx, col_idx)
            
            if task_id not in completed_tasks:
                task = (row_idx, col_idx, k_y_val, k_x_val, k_y_index, k_x_index, 
                        X_init, base_k, fi, t_span, t_eval)
                tasks.append(task)
                
    total_new_tasks = len(tasks)
    print(f"Total simulations remaining to run: {total_new_tasks}")
    
    if total_new_tasks > 0:
        # --- EXECUTE IN PARALLEL ---
        num_cores = max(1, mp.cpu_count() - 1) 
        print(f"Firing up {num_cores} CPU cores ")
        start_time = time.time()
        
        with mp.Pool(processes=num_cores) as pool:
            for i, result in enumerate(pool.imap_unordered(run_2d_simulation, tasks, chunksize=5)):
                row_idx, col_idx, alpha, r2 = result
                
                completed_tasks[(row_idx, col_idx)] = (alpha, r2)
                
                if (i + 1) % 50 == 0:
                    percent = ((i+1) / total_new_tasks) * 100
                    print(f"Completed {i + 1}/{total_new_tasks} remaining runs ({percent:.1f}%)...")
                    with open(checkpoint_file, "wb") as f:
                        pickle.dump(completed_tasks, f)
        
        with open(checkpoint_file, "wb") as f:
            pickle.dump(completed_tasks, f)
            
        elapsed = (time.time() - start_time) / 60
        print(f"\nFinished running new simulations in {elapsed:.2f} minutes!")
    else:
        print("All simulations were already completed! Proceeding to generate final grids...")

    # --- RECONSTRUCT FINAL GRIDS FROM THE CHECKPOINT ---
    print("\nCalculating final grids...")
    grid_shape = (len(k_change), len(k_change))
    alpha_grid = np.zeros(grid_shape)
    r2_grid = np.zeros(grid_shape)

    for (row, col), (alpha, r2) in completed_tasks.items():
        alpha_grid[row, col] = alpha
        r2_grid[row, col] = r2
    
    alpha_file = os.path.join(out_dir, f"alpha_{prefix}.pkl")
    r2_file = os.path.join(out_dir, f"r2_{prefix}.pkl")
    
    with open(alpha_file, "wb") as f:
        pickle.dump(alpha_grid, f)
        
    with open(r2_file, "wb") as f:
        pickle.dump(r2_grid, f)

    print(f"Saved to:\n - {alpha_file}\n - {r2_file}")

if __name__ == '__main__':
    mp.freeze_support()
    main()