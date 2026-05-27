from matplotlib import ticker
import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
from multiprocessing import Pool, cpu_count
import sys
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from src.fun_gilles import *

# --- Configuration ---
f = "examples/reactions_XYC.txt"
pkl_filename = "../data_archive/xyc_sweep_kar.pkl"
reps = 5
V0 = 100

initial_c_sweep = np.linspace(0, 500, 26) 
initial_c_sweep[0] = 1  # Avoid zero to prevent probability = 0
k_var = np.logspace(-4, 3, 8)
n_iterations = 1e5  



def run_single_condition(args):
    """Helper function to run the 5 reps for one specific condition."""
    c_val, k_val, f_path, v_init = args
    
    k_local = [1.0, 1.0, 1.0, k_val] 
    init_abundances = [c_val, 100, 100, 0, 0] 
    
    results = []
    for _ in range(reps):
        abundances, times, V = chemistry("Gillespie", n_iterations, f_path, 
                                        init_abundances, k_local, v_init)
        results.append(abundances[-1, -1] / V[-1])
    
    # NEW: Return BOTH the mean and the standard deviation
    return np.mean(results), np.std(results)

if __name__ == "__main__":
    
    if os.path.exists(pkl_filename):
        with open(pkl_filename, "rb") as pfile:
            master_data = pickle.load(pfile)
        print("Resuming from saved checkpoint...")
    else:
        master_data = {}

    for kv in k_var:
        k_str = f"{kv:.1e}"
        if k_str in master_data:
            print(f"Skipping k={k_str}, already calculated.")
            continue
        
        max_workers = 4 # Keep this at 4 to protect your RAM!
        print(f"Running sweep for k={k_str} using {max_workers} cores...")
        
        tasks = [(c, kv, f, V0) for c in initial_c_sweep]
        
        with Pool(processes=max_workers) as pool:
            # eq_values is now a list of tuples: [(mean1, std1), (mean2, std2), ...]
            eq_values = pool.map(run_single_condition, tasks)
        
        # NEW: Unpack the tuples into separate arrays
        means = np.array([val[0] for val in eq_values])
        stds = np.array([val[1] for val in eq_values])
        
        # Save as a dictionary inside the master_data
        master_data[k_str] = {'mean': means, 'std': stds}
        
        with open(pkl_filename, "wb") as pfile:
            pickle.dump(master_data, pfile)

    # --- Plotting ---
    plt.figure(figsize=(5,3))

    colors = plt.cm.rainbow(np.linspace(0, 1, len(k_var)))
    
    for i, kv in enumerate(k_var):
        k_label = f"{kv:.1e}"
        
        # Extract means and standard deviations
        means = master_data[k_label]['mean']
        stds = master_data[k_label]['std']
        
        # NEW: Use errorbar instead of plot
        plt.errorbar(initial_c_sweep/V0, means, yerr=stds, fmt='o-', 
                     color=colors[i], label=k_label, markersize=4, 
                     alpha=0.8, capsize=3, elinewidth=1.5)

    plt.xlabel("$[C]_0$", fontsize=14)
    plt.ylabel("$[XY]_{eq}$", fontsize=14)
    plt.tick_params(axis='both',  labelsize=13)
    plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    plt.legend(loc='center left',bbox_to_anchor=(1., 0.5), title="$k_{ar}$", fontsize=11, title_fontsize=12)
    plt.title("(A)", loc="left",x = -0.05, y=1, pad=10, fontsize=16)
    plt.tight_layout()
    plt.savefig("figures_TFM/gillespie_equilibrium_errors_kar.png", bbox_inches='tight', dpi=300)
    plt.show()