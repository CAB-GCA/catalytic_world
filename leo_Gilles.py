import numpy as np
import matplotlib.pyplot as plt
from fun_gilles import *

# Initialization:
file = "reactions_XY.txt" # M reactions
n_iterations = 1000
method = "Gillespie" # Gillespie or Deterministic
# Reaction constants:
k = [10,1,10,1] # len(k)= # de reacciones
# Volume:
V = 1
initial_food = 1000 # initial molecules number
food_molecules = 3

reactions = read_file(file)
species = obtain_species(reactions)


abundances, times = chemistry(method, n_iterations, reactions, food_molecules, 
                              initial_food, k, V)

# Representation
plt.figure()
colors = plt.cm.coolwarm(np.linspace(0, 1, len(species)))
plt.grid()

for i in range(len(species)):
    plt.plot(times, abundances[:, i], label=species[i], color=colors[i], alpha=0.9)

plt.grid(True, linestyle='--', alpha=0.3)
plt.xlabel("Time")
plt.ylabel("Abundances")
plt.legend()
plt.title(f"{method} time evolution: {n_iterations} iterations, k_a= {k[0]}")
    
plt.show()