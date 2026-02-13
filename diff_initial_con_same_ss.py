import matplotlib.pyplot as plt
import numpy as np
from fun_gilles import chemistry

# Setup parameters
file_path = "./examples/reactions_autocat.txt"
k_constants = [1]*8 # Example constants
total_time = 1e5
V_initial = 100

# Define different initial conditions
# IC1 and IC2 have the same total sum (100)
# IC3 has a different total sum (150)
initial_conditions = {
    "Group A (Sum=100) - Start 1": [80, 20], 
    "Group A (Sum=100) - Start 2": [20, 80],
    "Group B (Sum=150) - Start 3": [75, 75]
}
initial_condition = [100]*4 + [0]*4
plt.figure(figsize=(10, 6))

for label, ic in initial_conditions.items():
    initial_condition[4] = ic[0]
    initial_condition[5] = ic[1]

    # Using 'Deterministic' for clear visual convergence without stochastic noise
    abundances, times, V_nodes = chemistry(
        method='Protocell', 
        iterations=total_time, 
        file=file_path, 
        initial_food=initial_condition, 
        k=k_constants, 
        V=V_initial
    )
    
    # Calculate concentrations (X / V)
    # If V is constant in Deterministic, V_nodes is the scalar V
    concentrations = abundances.T / V_nodes 
    
    # Plot the first species for comparison
    plt.plot(times, concentrations[:, -1], label=f"{label} (Species 1)")

plt.title("Convergence to Steady State based on Total Concentration")
plt.xlabel("Time")
plt.ylabel("Concentration [X]")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

