import numpy as np
import matplotlib.pyplot as plt
from fun_gilles import *

# Initialization:
file = "reactions_XYC.txt" # M reactions
n_iterations = 10000
method = "Gillespie"
# Reaction constants:
k = [0.01]*4
# Volume:
V = 100
initial_food = 50 # Initial population number
food_molecules = 3

reactions = read_file(file)
species = obtain_species(reactions)

abundances, times = chemistry(method, n_iterations, reactions, food_molecules, initial_food, k, V)

# Representation
plt.figure()
colors = plt.cm.coolwarm(np.linspace(0, 1, len(species)))
plt.grid()

for i in range(len(species)):
    plt.plot(times, abundances[:, i], label=species[i], color=colors[i], alpha=0.8)

plt.legend()
    
plt.show()