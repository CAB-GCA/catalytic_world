# -*- coding: utf-8 -*-
"""
Created on Thu May 21 11:41:56 2026

@author: cvzad
"""

import os
import sys
import pickle

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from src.fun_gilles import chemistry

def save_single_run_to_pickle(ab, t, v, filename):
    """
    Saves the simulation tuple into a list to match the format 
    expected by load_list_pickle() in internship_figures.py
    """
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "wb") as file:
        pickle.dump([(ab, t, v)], file)
    print(f"Successfully generated and saved: {filename}")


if __name__ == "__main__":
    print("Starting generation of Internship Report .pkl files...\n")
    
    # Define the output directory
    out_dir = os.path.join(project_root, "practicas")
    
    # Common Parameters
    V_init = 1000.0  # Assumed standard volume scale to convert concentrations to abundances

    # ============================================================
    # 1. XYC SYSTEMS (Figure 4)
    # ============================================================
    print("Generating XYC simulations...")
    xyc_file = os.path.join(project_root, "examples", "reactions_XYC_food.txt")
    
    # Fig 4a: Constant Food (Protocell)

    xyc_proto_ab0 = [0.25 * V_init, 0.25 * V_init, 1.0 * V_init, 0, 0]
    
    ab, t, v = chemistry(
        method="Protocell", iterations=20000, file=xyc_file,
        initial_food=xyc_proto_ab0, k=[1.0]*4, V=V_init
    )
    save_single_run_to_pickle(ab, t, v, os.path.join(out_dir, "XYC_protocell.pkl"))
    save_single_run_to_pickle(ab, t, v, os.path.join(out_dir, "xyc_gillespie.pkl")) 

    # Deterministic counterpart for 2B
    ab_det, t_det, v_det = chemistry(
        method="Deterministic", iterations=100, file=xyc_file,
        initial_food=xyc_proto_ab0, k=[1.0]*6, V=V_init
    )
    save_single_run_to_pickle(ab_det.T, t_det, v_det, os.path.join(out_dir, "xyc_deterministic.pkl"))


    # ============================================================
    # 2. ABCD (AUTOCAT) SYSTEMS - GOOD RUNS (Figure 5a & 5b)
    # ============================================================
    print("\nGenerating ABCD 'Good' simulations (ab0 = 1)...")
    autocat_file = os.path.join(project_root, "examples", "reactions_autocat.txt")
    
    # Fig 5a: Constant Food, Low Catalyst
    abcd_good_ab0 = [1.0 * V_init]*4 + [1.0 * V_init] + [0]*3
    
    ab, t, v = chemistry(
        method="Protocell", iterations=100000, file=autocat_file,
        initial_food=abcd_good_ab0, k=[1.0]*8, V=V_init, threshold=0
    )
    save_single_run_to_pickle(ab, t, v, os.path.join(out_dir, "abcd_protocell.pkl"))

    # Fig 5b: Constant Flux, Low Catalyst
    abcd_flux_good_ab0 = [0]*4 + [1.0 * V_init] + [0]*3
    k_flux_good = [1.0]*8 + [0.8]*4 
    
    ab, t, v = chemistry(
        method="Gillespie", iterations=100000, file=autocat_file,
        initial_food=abcd_flux_good_ab0, k=k_flux_good, V=V_init
    )
    save_single_run_to_pickle(ab, t, v, os.path.join(out_dir, "abcd_ctflux.pkl"))


    # ============================================================
    # 3. ABCD (AUTOCAT) SYSTEMS - BAD RUNS (From Figure 5c & 5d)
    # ============================================================
    print("\nGenerating ABCD 'Bad' simulations (ab0 = 6)...")
    
    # Fig 5c: Constant Food, High Catalyst
    abcd_bad_ab0 = [1.0 * V_init]*4 + [6.0 * V_init] + [0]*3
    
    ab, t, v = chemistry(
        method="Protocell", iterations=100000, file=autocat_file,
        initial_food=abcd_bad_ab0, k=[1.0]*8, V=V_init, threshold=0
    )
    save_single_run_to_pickle(ab, t, v, os.path.join(out_dir, "abcd_protocell_bad.pkl"))

    # Fig 5d: Constant Flux, High Catalyst
    abcd_flux_bad_ab0 = [0]*4 + [6.0 * V_init] + [0]*3
    
    ab, t, v = chemistry(
        method="Gillespie", iterations=500000, file=autocat_file,
        initial_food=abcd_flux_bad_ab0, k=k_flux_good, V=V_init
    )
    save_single_run_to_pickle(ab, t, v, os.path.join(out_dir, "abcd_ctflux_bad.pkl"))

    print("\nAll required .pkl files have been generated in the 'practicas/' folder!")