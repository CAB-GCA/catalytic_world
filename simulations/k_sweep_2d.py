# TO GENERATE THE PKL FILES RUN "k_sweep_2d_pkl_generator.py"


import numpy as np
import os
import sys
import pickle


current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from src.plotting import (set_publication_style, 
                          plot_precalculated_alpha, 
                          plot_r2_map, 
                          plot_both_2d_sweeps)


def load_alpha_grid(filepath):
    """Safely loads a pre-calculated 2D grid pickle."""
    try:
        with open(filepath, "rb") as f:
            return pickle.load(f)
    except FileNotFoundError:
        print(f"Error: Could not find {filepath}")
        return np.array([])


if __name__ == "__main__":
    
    set_publication_style()
    print("Generating 2D Parameter Sweep Figures...\n")
    
    # Common logspace used for both sets
    k_change = np.logspace(-5, 7, 25)

    # ============================================================
    # 1. SWEEP: k2 vs k4 (using files alpha_sweep_k1_k3 / r2_sweep_k1_k3)
    # ============================================================
    print("Processing k2 vs k4 sweep...")
    path_alpha_1_3 = os.path.join(project_root, "data_archive", "alpha_sweep_k1_k3_deterministic.pkl")
    path_r2_1_3    = os.path.join(project_root, "data_archive", "r2_sweep_k1_k3_deterministic.pkl")
    
    k2_vs_k4 = load_alpha_grid(path_alpha_1_3)
    r2_grid_1_3 = load_alpha_grid(path_r2_1_3)
    
    if k2_vs_k4.size > 0 and r2_grid_1_3.size > 0:
        plot_precalculated_alpha(k2_vs_k4, k_change, k_in=4, k_out=2, figname="alpha_k2_vs_k4")
        plot_r2_map(r2_grid_1_3, k_change, 4, 2, figname="r2_k2_vs_k4")
        plot_both_2d_sweeps(k2_vs_k4, r2_grid_1_3, k_change, k_in=4, k_out=2, figname="alpha_and_r2_k2_vs_k4")

    # ============================================================
    # 2. SWEEP: k1 vs k3 (using files alpha_sweep_k0_k2 / r2_sweep_k0_k2)
    # ============================================================
    print("\nProcessing k1 vs k3 sweep...")
    path_alpha_0_2 = os.path.join(project_root, "data_archive", "alpha_sweep_k0_k2.pkl")
    path_r2_0_2    = os.path.join(project_root, "data_archive", "r2_sweep_k0_k2.pkl")
    
    k1_vs_k3 = load_alpha_grid(path_alpha_0_2)
    r2_grid_0_2 = load_alpha_grid(path_r2_0_2)
    
    if k1_vs_k3.size > 0 and r2_grid_0_2.size > 0:
        plot_precalculated_alpha(k1_vs_k3, k_change, k_in=3, k_out=1, figname="alpha_k1_vs_k3")
        plot_r2_map(r2_grid_0_2, k_change, 3, 1, figname="r2_k1_vs_k3")
        plot_both_2d_sweeps(k1_vs_k3, r2_grid_0_2, k_change, k_in=3, k_out=1, figname="alpha_and_r2_k1_vs_k3")
        
    print("\nDone! All figures saved in 'figures_TFM/'")