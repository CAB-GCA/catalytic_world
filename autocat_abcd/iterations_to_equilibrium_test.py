import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# --- Path setup ---
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from fun_gilles import *

f = "examples/reactions_autocat.txt"
n_iterations = int(1e10) 

initial_food = [10, 10, 10, 10] + [2, 0] + [0, 0]
v0 = 10
k = [1, 0] * 4 + [0] * 4 # nutrient influx is set to zero and reactions are irreversible. the system will stop when all nutrients are consumed
volume_factor = np.linspace(0, 1000, 11)
volume_factor[0] = 1 

# --- Data Containers ---
volumes = []
iterations_count = []

for i in volume_factor:
    current_volume = v0 * i
    print(f"Running Gillespie simulation with volume: {current_volume}")
    
    initial_food_scaled = [x * i for x in initial_food]
    
    # Run simulation
    a, t, v = chemistry("Gillespie", n_iterations, f, initial_food_scaled, k, current_volume)    
    
    # Store results
    volumes.append(current_volume)
    iterations_count.append(len(t))
    
    print(f"Number of reactions occurred: {len(t)}")

# --- Plotting ---
plt.figure(figsize=(5,4))
plt.plot(volumes, iterations_count, marker='o', linestyle='-', color='teal', linewidth=2)

plt.title("Computational Cost vs. System Size", fontsize=14)
plt.xlabel("Volume ($V$)", fontsize=12)
plt.ylabel("Number of Iterations (Reactions)", fontsize=12)
plt.grid(True, which="both", ls="--", alpha=0.5)

# plt.xscale('log')
# plt.yscale('log')

plt.tight_layout()
plt.show()


x = np.array(volumes)
y = np.array(iterations_count)

# 2. Perform Linear Regression
slope, intercept, r_value, p_value, std_err = linregress(x, y)

# 3. Define a function to predict future iterations
def predict_iterations(target_volume):
    return int(slope * target_volume + intercept)

print(f"Regression Results:")
print(f"  Slope: {slope:.4f}")
print(f"  Intercept: {intercept:.4f}")
print(f"  R-squared: {r_value**2:.4f}")

# Example prediction
example_v = 5000
predicted_n = predict_iterations(example_v)
print(f"\nPredicted iterations for Volume {example_v}: {predicted_n}")

# 4. Visualization
plt.figure(figsize=(8, 5))
plt.scatter(x, y, color='teal', label='Simulated Data')
plt.plot(x, slope*x + intercept, color='crimson', label=f'Fit: y={slope:.2f}x + {intercept:.2f}')

plt.xlabel('Volume (V)')
plt.ylabel('Total Iterations')
plt.title('Predicting Simulation Cost via Linear Regression')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()