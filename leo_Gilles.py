import numpy as np
import matplotlib.pyplot as plt
from fun_gilles import *

# Initialization:
file = "reactions_L_4_d_1.txt" # M reactions
n_iterations = 10
method = "Deterministic"
# Reaction constants:
k = [0.1]*12
# Volume:
V = 10
initial_food = 500 # Initial population number
food_molecules = 3

reactions = read_file(file)
species = obtain_species(reactions)

abundances, times = chemistry(method, n_iterations, reactions, food_molecules, initial_food, k, V)

# Representation
plt.figure()
colors = plt.cm.Spectral(np.linspace(0, 1, len(species)))
plt.grid()

for i in range(len(species)):
    plt.plot(times, abundances[:, i], label=species[i], color=colors[i], alpha=0.8)

plt.legend()
    
plt.show()