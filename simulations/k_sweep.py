# TO GENERATE THE PKL FILES RUN "k_sweep_pkl_generator.py"
import os
import sys
import pickle

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from src.plotting import set_publication_style, plot_combined_k_sweep_and_volumes

def load_streamed_pickle(filename):
    """Robust function to load streamed dictionaries from a pickle file."""
    data = {}
    try:
        with open(filename, "rb") as file:
            while True:
                try:
                    chunk = pickle.load(file)
                    data.update(chunk)
                except EOFError:
                    break
        print(f"Successfully loaded {len(data)} simulation results from {os.path.basename(filename)}.")
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    return data

if __name__ == "__main__":
    
    set_publication_style()
    print("Generating 1D Kinetic Parameter Sweep Figures...\n")

    # ============================================================
    # 1. SWEEP FOR k2 (i=1)
    # ============================================================
    i_k2 = 1
    file_k2 = os.path.join(project_root, "data_archive", f"barrido_k{i_k2}.pkl")
    print(f"Processing sweep for k_{i_k2 + 1}...")
    
    results_k2 = load_streamed_pickle(file_k2)
    if results_k2:
        plot_combined_k_sweep_and_volumes(
            results_k2, 
            k_i=i_k2, 
            xlim=50, 
            lim=(5e-2, 1e4), 
            figname=f"alpha_and_vol_k_{i_k2 + 1}"
        )

    # ============================================================
    # 2. SWEEP FOR k1 (i=0)
    # ============================================================
    i_k1 = 0
    file_k1 = os.path.join(project_root, "data_archive", f"barrido_k{i_k1}.pkl")
    print(f"\nProcessing sweep for k_{i_k1 + 1}...")
    
    results_k1 = load_streamed_pickle(file_k1)
    if results_k1:
        plot_combined_k_sweep_and_volumes(
            results_k1, 
            k_i=i_k1, 
            xlim=120, 
            lim=(5e-4, 1e2), 
            figname=f"alpha_and_vol_k_{i_k1 + 1}"
        )
        
    print("\nDone! All figures saved in 'figures_TFM/'")