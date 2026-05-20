from scipy.integrate import solve_ivp
import numpy as np
import pickle


k_values= [1e13,1e16,1e-1,1,1e13,1e16,1e-1,1]

def odes_intermediates(t, X):
    V, ab, cd, acd, cab = X
    a, b, c, d = 1, 1, 1, 1
    k1, k2, k3, k4, k5, k6, k7, k8 = k_values
    fi = 0.2
    dVdt = (V/fi) * (k3 * b * acd + k7 * d * cab - k4 * ab * cd - k8 * ab * cd)
    dabdt = k3 * b * acd + k6 * cab + k7 * cab * d - k4 * ab * cd - k5 * ab * c - k8 * cd * ab - (ab / V) * dVdt
    dcddt = k3 * b * acd + k7 * cab * d + k2 * acd - k4 * ab * cd - k1 * a * cd - k8 * cd * ab - (cd / V) * dVdt
    dacddt = k4 * ab * cd + k1 * a * cd - k3 * b * acd - k2 * acd - (acd / V) * dVdt
    dcabdt = k5 * ab * c + k8 * cd * ab - k6 * cab - k7 * cab * d - (cab / V) * dVdt
    return [dVdt, dabdt, dcddt, dacddt, dcabdt]
    

t_span = (0, 1000)
y0 = [1000, 0.2, 0, 0, 0]
sol = solve_ivp(odes_intermediates, t_span, y0, t_eval=np.linspace(0, 1000, 500), method='BDF')

print("ODE solution computed successfully.")
with open("third_order/odes_intermediates_solution.pkl", "ab") as f:
    pickle.dump({k_values[2] : sol}, f)

