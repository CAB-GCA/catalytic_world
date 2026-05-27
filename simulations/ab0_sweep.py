###############################################################################################
# FOR AB0_SWEEP THE WORKFLOW IS THE FOLLOWING:
    # 1. GENERATE THE PKL FILES USING ab0_sweep_pkl_generator.py SETTING THE DESIRED PARAMETERS
    # 2. RUN ab0_sweep.py TO GET THE FIGURES AND STATS
    # 3. TO CHANGE THE TITLES OF THE FIGURE (WHICH MAY NOT CORRESPOND TO THE CHOSEN 
    #                                        PARAMETERS) GO TO ../src/plotting.py
###############################################################################################
import numpy as np
import pandas as pd
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)
from src.fun_gilles import *
from src.plotting import plot_combined_AB, set_publication_style, load_streamed_pickle, extract_stats_from_replicates
import pickle


if __name__ == "__main__":
    # Setup styles
    set_publication_style()
    
    # Load data for the triple sweep
    print("Loading data...")
    results_k1_low = load_streamed_pickle("../data_archive/barrido_ab0_k1_0.0001_parallel.pkl")
    results_k1_med = load_streamed_pickle("../data_archive/barrido_ab0_k1_1_local.pkl")
    results_k1_high = load_streamed_pickle("../data_archive/barrido_ab0_k1_10000.0_parallel.pkl")
    
    # Optional: If you want to print the DataFrames to check the exact numbers before plotting
    print("\nStats for k1 = 10^-4:")
    print(extract_stats_from_replicates(results_k1_low))
    
    # Generate the massive 2x3 combined figure
    print("\nGenerating Figure...")
    plot_combined_AB(
        results_list=[results_k1_low, results_k1_med, results_k1_high],
        xlim=[18, 50, 60], 
        fig_name="ab0_combined_final"
    )
    print("Done!")