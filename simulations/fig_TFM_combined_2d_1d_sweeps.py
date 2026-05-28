import numpy as np
import os
import sys
import pickle

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from src.plotting import set_publication_style, plot_master_4panel_figure

def load_alpha_grid(filepath):
    """Safely loads a pre-calculated 2D grid pickle."""
    try:
        with open(filepath, "rb") as f:
            return pickle.load(f)
    except FileNotFoundError:
        print(f"Error: Could not find {os.path.basename(filepath)}")
        return np.array([])

def load_streamed_pickle(filepath):
    """Robust function to load streamed dictionaries from a 1D sweep pickle file."""
    data = {}
    try:
        with open(filepath, "rb") as file:
            while True:
                try:
                    chunk = pickle.load(file)
                    data.update(chunk)
                except EOFError:
                    break
        return data
    except FileNotFoundError:
        print(f"Error: Could not find {os.path.basename(filepath)}")
        return {}

if __name__ == "__main__":
    set_publication_style()
    print("Generating Master 4-Panel Figures...\n")
    
    k_change_2d = np.logspace(-5, 7, 25)

    # ============================================================
    # MASTER FIGURE 1: Focus on k1 (k1 vs k3 2D map + k1 1D sweep)
    # ============================================================
    print("Building Master Figure for k1...")
    path_alpha_0_2 = os.path.join(project_root, "data_archive", "alpha_sweep_k0_k2.pkl")
    path_r2_0_2    = os.path.join(project_root, "data_archive", "r2_sweep_k0_k2.pkl")
    path_k1_1d     = os.path.join(project_root, "data_archive", "barrido_k0.pkl")
    
    k1_vs_k3 = load_alpha_grid(path_alpha_0_2)
    r2_grid_0_2 = load_alpha_grid(path_r2_0_2)
    results_k1 = load_streamed_pickle(path_k1_1d)
    
    if k1_vs_k3.size > 0 and results_k1:
        plot_master_4panel_figure(
            alpha_grid=k1_vs_k3, r2_grid=r2_grid_0_2, k_values_2d=k_change_2d, k_in=3, k_out=1,
            results_1d=results_k1, k_i=0, 
            xlim_1d=120, lim_1d=(5e-4, 1e2),
            figname="master_sweep_k1"
        )

    # ============================================================
    # MASTER FIGURE 2: Focus on k2 (k2 vs k4 2D map + k2 1D sweep)
    # ============================================================
    print("Building Master Figure for k2...")
    path_alpha_1_3 = os.path.join(project_root, "data_archive", "alpha_sweep_k1_k3.pkl")
    path_r2_1_3    = os.path.join(project_root, "data_archive", "r2_sweep_k1_k3.pkl")
    path_k2_1d     = os.path.join(project_root, "data_archive", "barrido_k1.pkl")
    
    k2_vs_k4 = load_alpha_grid(path_alpha_1_3)
    r2_grid_1_3 = load_alpha_grid(path_r2_1_3)
    results_k2 = load_streamed_pickle(path_k2_1d)
    
    if k2_vs_k4.size > 0 and results_k2:
        plot_master_4panel_figure(
            alpha_grid=k2_vs_k4, r2_grid=r2_grid_1_3, k_values_2d=k_change_2d, k_in=4, k_out=2,
            results_1d=results_k2, k_i=1, 
            xlim_1d=50, lim_1d=(5e-2, 1e4),
            figname="master_sweep_k2"
        )
        
    print("Done! Master figures saved to 'figures_TFM/'")

