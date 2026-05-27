# ============================================================
# Refactored: internship_figures.py
# ============================================================
import numpy as np
import sys
import os
import pickle

# --- ROBUST IMPORT FIX ---
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from src.fun_gilles import chemistry, read_file, obtain_species
from src.plotting import (set_publication_style, 
                          plot_internship_abundances, 
                          plot_internship_conc_and_vol)

def load_list_pickle(filename):
    data = []
    try:
        with open(filename, "rb") as file:
            while True:
                try:
                    chunk = pickle.load(file)
                    # If the chunk is already a list, extend the data. 
                    # If it's just a tuple, append it.
                    if isinstance(chunk, list):
                        data.extend(chunk)
                    else:
                        data.append(chunk)
                except EOFError:
                    break
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    
    return data

if __name__ == "__main__":
    set_publication_style()
    print("Generating Internship Report Figures...\n")

    # ============================================================
    # 1. XY SYSTEM (Direct Calculation)
    # ============================================================
    print("Plotting XY System...")
    file_XY = os.path.join(project_root, "examples", "reactions_XY.txt")
    species_XY = obtain_species(read_file(file_XY))
    
    ab_xy, t_xy, v_xy = chemistry(
        method="Gillespie", iterations=1200, file=file_XY,
        initial_food=[1000, 1000, 0], k=[1, 1], V=1
    )
    plot_internship_abundances(t_xy, ab_xy, species_XY, fig_name="XY_Gillespie")


    # ============================================================
    # 2. XYC SYSTEM (Pre-calculated Pickles)
    # ============================================================
    print("Plotting XYC Systems...")
    xyc_food_file = os.path.join(project_root, "examples", "reactions_XYC_food.txt")
    species_xyc = obtain_species(read_file(xyc_food_file))
    
    # 2A. XYC Gillespie (Cut off last 200 to match deterministic length)
    xyc_gill_path = os.path.join(project_root, "pkl_for_internship", "xyc_gillespie.pkl")
    if os.path.exists(xyc_gill_path):
        ab, t, v = load_list_pickle(xyc_gill_path)[0]
        plot_internship_abundances(t[:-200], ab[:-200, :], species_xyc, fig_name="XYC_Gillespie")

    # 2B. XYC Deterministic
    xyc_det_path = os.path.join(project_root, "pkl_for_internship", "xyc_deterministic.pkl")
    if os.path.exists(xyc_det_path):
        ab, t, v = load_list_pickle(xyc_det_path)[0]
        plot_internship_abundances(t, ab.T, species_xyc, fig_name="XYC_Deterministic")

    # 2C. XYC Protocell
    xyc_proto_path = os.path.join(project_root, "pkl_for_internship", "XYC_protocell.pkl")
    if os.path.exists(xyc_proto_path):
        ab, t, v = load_list_pickle(xyc_proto_path)[0]
        plot_internship_conc_and_vol(t, ab, v, species_xyc, fig_name="XYC_Protocell")


    # ============================================================
    # 3. ABCD (AUTOCAT) SYSTEM - GOOD RUNS
    # ============================================================
    print("Plotting ABCD Good Systems...")
    autocat_file = os.path.join(project_root, "examples", "reactions_autocat.txt")
    species_abcd = obtain_species(read_file(autocat_file))
    
    # 3A. ABCD Protocell
    abcd_proto_path = os.path.join(project_root, "pkl_for_internship", "abcd_protocell.pkl")
    if os.path.exists(abcd_proto_path):
        ab, t, v = load_list_pickle(abcd_proto_path)[0]
        plot_internship_conc_and_vol(t, ab, v, species_abcd, fig_name="ABCD_Protocell_Good")

    # 3B. ABCD Constant Flux
    abcd_ctflux_path = os.path.join(project_root, "pkl_for_internship", "abcd_ctflux.pkl")
    if os.path.exists(abcd_ctflux_path):
        ab, t, v = load_list_pickle(abcd_ctflux_path)[-1]
        plot_internship_conc_and_vol(t, ab, v, species_abcd, fig_name="ABCD_ConstantFlux_Good")


    # ============================================================
    # 4. ABCD (AUTOCAT) SYSTEM - BAD RUNS
    # ============================================================
    print("Plotting ABCD Bad Systems...")
    
    # 4A. ABCD Protocell Bad
    abcd_proto_bad_path = os.path.join(project_root, "pkl_for_internship", "abcd_protocell_bad.pkl")
    if os.path.exists(abcd_proto_bad_path):
        ab, t, v = load_list_pickle(abcd_proto_bad_path)[-1]
        plot_internship_conc_and_vol(t, ab, v, species_abcd, fig_name="ABCD_Protocell_Bad")

    # 4B. ABCD Constant Flux Bad
    abcd_ctflux_bad_path = os.path.join(project_root, "pkl_for_internship", "abcd_ctflux_bad.pkl")
    if os.path.exists(abcd_ctflux_bad_path):
        ab, t, v = load_list_pickle(abcd_ctflux_bad_path)[-1]
        plot_internship_conc_and_vol(t, ab, v, species_abcd, fig_name="ABCD_ConstantFlux_Bad")
        
    print("\nDone! All figures saved in 'figures_internship/'")