import numpy as np
import matplotlib.pyplot as plt
from fun_gilles import *

# Initialization:
file = "reactions_XYC.txt" # M reactions
n_iterations = 210
method = "Gillespie" # Gillespie or Deterministic
# Reaction constants:
k = [1,.01]*2 # len(k)= # de reacciones
# Volume:
V = 0.1
initial_food = [120,100,120] # initial molecules number
food_molecules = 3

reactions = read_file(file)
species = obtain_species(reactions)


abundances, times = chemistry(method, n_iterations, reactions, food_molecules, 
                              initial_food, k, V)

it= n_iterations
# Representation
plt.figure()
colors = plt.cm.coolwarm(np.linspace(0, 1, len(species)))
plt.grid()

for i in range(len(species)):
    plt.plot(times[:it], abundances[:it, i], label=species[i], color=colors[i], alpha=0.8)

plt.grid(True, linestyle='--', alpha=0.3)
plt.xlabel("Time")
plt.ylabel("Abundances")
plt.legend()
plt.title(f"{method} time evolution. V={V}")

    
plt.show()