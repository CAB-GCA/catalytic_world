
import numpy as np
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from src.fun_gilles import *
from src.plotting import (set_publication_style, 
                          plot_concentration_evolution, 
                          plot_comparison_concentration, 
                          plot_comparison_conc_and_vol)

if __name__ == "__main__":
    
    set_publication_style()
    print("Executing Time-Evolution Simulations...")

    # ============================================================
    # 1. SIMPLE XYC SYSTEM
    # ============================================================
    print("Running XYC model...")
    xyc_file = os.path.join(project_root, "examples", "reactions_XYC.txt")
    xyc_species = ["c", "x", "y", "xc", "xy"]
    xyc_species = [word.upper() for word in xyc_species]

    a_xyc, t_xyc, v_xyc = chemistry(
        method="Gillespie", iterations=5e4, file=xyc_file,
        initial_food=[100, 1000, 1000, 0, 0], k=[1]*4, V=1000
    )
    
    a_det, t_det, v_det = chemistry(
        method="Deterministic", iterations=300, file=xyc_file,
        initial_food=[100, 1000, 1000, 0, 0], k=[1]*4, V=1000
    )

    plot_concentration_evolution(t_xyc, a_xyc, v_xyc, xyc_species, fig_name="XYC_time_evolution")
    plot_comparison_concentration((a_det.T/v_det, t_det), (a_xyc, t_xyc, v_xyc), 
                                  xyc_species, fig_name="XYC_time_evolution_comparison", xlim=300)

    # ============================================================
    # 2. AUTOCATALYTIC PROTOCELL SYSTEM 
    # ============================================================
    print("Running Autocatalytic model (Stochastic)...")
    autocat_file = os.path.join(project_root, "examples", "reactions_autocat.txt")
    species = obtain_species(read_file(autocat_file))
    species = [word.upper() for word in species]

    
    a, t, v = chemistry(
        method="Protocell", file=autocat_file, iterations=2e5,
        initial_food=[1000]*4 + [1000]*1 + [0]*3, k=[1]*8, V=1000, threshold=0
    )
    print(f"Final Stochastic Time: {t[-1]}")
    plot_concentration_evolution(t, a, v, species, fig_name="autocat_time_evolution")

    # ============================================================
    # 3. AUTOCATALYTIC ODE INTEGRATION 
    # ============================================================
    print("Running Autocatalytic model (Deterministic)...")
    def odes_intermediates(t, X, K_intermediates):
        V, ab, cd, acd, cab = X
        a, b, c, d = 1, 1, 1, 1
        k1, k2, k3, k4, k5, k6, k7, k8 = K_intermediates
        fi = 1
        dVdt = (V/fi) * (k3 * b * acd + k7 * d * cab - k4 * ab * cd - k8 * ab * cd)
        dabdt = k3 * b * acd + k6 * cab + k7 * cab * d - k4 * ab * cd - k5 * ab * c - k8 * cd * ab - (ab / V) * dVdt
        dcddt = k3 * b * acd + k7 * cab * d + k2 * acd - k4 * ab * cd - k1 * a * cd - k8 * cd * ab - (cd / V) * dVdt
        dacddt = k4 * ab * cd + k1 * a * cd - k3 * b * acd - k2 * acd - (acd / V) * dVdt
        dcabdt = k5 * ab * c + k8 * cd * ab - k6 * cab - k7 * cab * d - (cab / V) * dVdt
        return [dVdt, dabdt, dcddt, dacddt, dcabdt]
        
    t_end = 17.4
    t_span = (0, t_end)
    y0 = [1000, 1, 0, 0, 0]
    K_intermediates = [1]*8
    
    sol = solve_ivp(odes_intermediates, t_span, y0, t_eval=np.linspace(0, t_end, 5000), 
                    method='BDF', args=(K_intermediates,))
    
    # Plot final comparison
    print("Generating final comparison figure...")
    plot_comparison_conc_and_vol(
        (sol.y[1:], sol.t, sol.y[0]), 
        (a[:, 4:], t, v), 
        species[4:], 
        fig_name="Autocatalytic_time_evolution_comparison"
    )
    print("All figures successfully generated and saved!")