from scipy.integrate import solve_ivp
import numpy as np
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from src.plotting import (set_publication_style, 
                          plot_deterministic_dual_comparison, 
                          plot_volume_sweep_dual_comparison)

# ============================================================
# ODE DEFINITIONS
# ============================================================

def odes_thirdorder(t, X, K_thirdorder):
    V, ab, cd = X
    a, b, c, d = 1, 1, 1, 1
    k1, k2, k3, k4 = K_thirdorder
    fi = 0.2
    
    dVdt = (V/fi) * (k1 * a * b * cd - k2 * ab * cd + k3 * ab * c * d - k4 * ab * cd)
    dabdt = k1 * a * b * cd - k2 * ab * cd - (ab / V) * dVdt 
    dcddt = k3 * ab * c * d - k4 * ab * cd - (cd / V) * dVdt
    return [dVdt, dabdt, dcddt]

def odes_intermediates(t, X, K_intermediates):
    V, ab, cd, acd, cab = X
    a, b, c, d = 1, 1, 1, 1
    k1, k2, k3, k4, k5, k6, k7, k8 = K_intermediates
    fi = 0.2
    
    dVdt = (V/fi) * (k3 * b * acd + k7 * d * cab - k4 * ab * cd - k8 * ab * cd)
    dabdt = k3 * b * acd + k6 * cab + k7 * cab * d - k4 * ab * cd - k5 * ab * c - k8 * cd * ab - (ab / V) * dVdt
    dcddt = k3 * b * acd + k7 * cab * d + k2 * acd - k4 * ab * cd - k1 * a * cd - k8 * cd * ab - (cd / V) * dVdt
    dacddt = k4 * ab * cd + k1 * a * cd - k3 * b * acd - k2 * acd - (acd / V) * dVdt
    dcabdt = k5 * ab * c + k8 * cd * ab - k6 * cab - k7 * cab * d - (cab / V) * dVdt
    return [dVdt, dabdt, dcddt, dacddt, dcabdt]


if __name__ == "__main__":
    set_publication_style()
    print("Running Kinetic Equivalence Analysis...\n")
    
    t_end = 12
    t_span = (0, t_end)
    t_eval = np.linspace(0, t_end, 500)
    
    # ============================================================
    # 1. DIRECT COMPARISON (DUAL GRID)
    # ============================================================
    print("Calculating direct comparison (Equivalent vs Broken)...")
    
    # Third-Order Base
    y0_to = [1000, 0.2, 0]
    K_to_base = [1, 1, 1, 1]
    sol_to = solve_ivp(odes_thirdorder, t_span, y0_to, t_eval=t_eval, args=(K_to_base,))
    
    y0_int = [1000, 0.2, 0, 0, 0]
    
    # Intermediates (EQUIVALENT - Rapid Equilibrium)
    K_int_eq = [1e4, 1e7, 1e3, 1, 1e4, 1e7, 1e3, 1]
    sol_int_eq = solve_ivp(odes_intermediates, t_span, y0_int, t_eval=t_eval, method='BDF', args=(K_int_eq,))
    
    # Intermediates (BROKEN - Slow intermediate formation, partial equilibrium fails)
    K_int_br = [1]*8
    sol_int_br = solve_ivp(odes_intermediates, t_span, y0_int, t_eval=t_eval, method='BDF', args=(K_int_br,))
    
    plot_deterministic_dual_comparison(
        sol_to.t, sol_to.y[1:3], sol_to.y[0], ['AB', 'CD'],
        sol_int_eq.t, sol_int_eq.y[1:], sol_int_eq.y[0],
        sol_int_br.t, sol_int_br.y[1:], sol_int_br.y[0], ['AB', 'CD', 'ACD', 'CAB'],
        figname="dual_equivalence_comparison"
    )

    # ============================================================
    # 2. PARAMETER SWEEP COMPARISON (DUAL GRID)
    # ============================================================
    print("Calculating parameter sweep (Equivalent vs Broken)...")
    
    t_end = 4
    t_span = (0, t_end)
    t_eval = np.linspace(0, t_end, 500)
    k_values = np.logspace(-3, 1, 9)
    results_to, results_int_eq, results_int_br = [], [], []
    
    for i in k_values:
        # Third order (Sweep k1)
        K_to_sweep = [i, 1, 1, 1]
        sol2 = solve_ivp(odes_thirdorder, t_span, y0_to, t_eval=t_eval, method='BDF', args=(K_to_sweep,))
        results_to.append((sol2.t, sol2.y[0])) 
        
        # Intermediates (EQUIVALENT)
        K_int_eq_sweep = [1e4, 1e7, 1, 1, 1e4, 1e7, 1e3, 1]
        K_int_eq_sweep[2] = i * K_int_eq_sweep[1] / K_int_eq_sweep[0]
        sol_eq = solve_ivp(odes_intermediates, t_span, y0_int, t_eval=t_eval, method='BDF', args=(K_int_eq_sweep,))
        results_int_eq.append((sol_eq.t, sol_eq.y[0]))
        
        # Intermediates (BROKEN)
        # Apply the exact same derivation formula, but the physical assumptions are violated by the low base rates.
        K_int_br_sweep = [1]*8
        K_int_br_sweep[2] = i * K_int_br_sweep[1] / K_int_br_sweep[0] 
        sol_br = solve_ivp(odes_intermediates, t_span, y0_int, t_eval=t_eval, method='BDF', args=(K_int_br_sweep,))
        results_int_br.append((sol_br.t, sol_br.y[0]))
        
    plot_volume_sweep_dual_comparison(results_to, results_int_eq, results_int_br, k_values, figname="dual_sweep_comparison")
    print("Done! Figures saved.")