import numpy as np
import matplotlib.pyplot as plt
from fun_gilles import *

# Initialization:
file = "reactions_XYC.txt" # M reactions
method = "Deterministic" # Gillespie or Deterministic

# Reaction constants:
k = [0,10,100,10] # len(k)= # de reacciones
# Volume:
V = 1

# condiciones iniciales
initial_food = [0,10,20] # initial molecules number
food_molecules = 3

# obtener reacciones y especies:
reactions = read_file(file)
species = obtain_species(reactions)

# Different values for initial C concentration and k_a:
initial_c = np.round(np.linspace(1,25,25))
k_a = np.logspace(1, 5, 5)
# To make sure the system reaches a "steady state" the number of iterations must 
# be different for each k_a
n_iterations = [20000,2500,820,600,600,600]

mean_equilibrium = np.zeros((len(initial_c)))
colors = plt.cm.Spectral(np.linspace(0, 1, len(k_a)))

for j in range(len(k_a)):
    k[0] = k_a[j]
    
    for i in range(len(initial_c)):
        initial_food[0] = initial_c[i]
        abundances, times = chemistry(method, n_iterations[j], reactions, 
                                      food_molecules, initial_food, k, V)
        mean_equilibrium[i] = np.mean(abundances[-min(50, len(abundances)):, -1])
    
    plt.plot(initial_c, mean_equilibrium, label=f"k_a={k_a[j]:.1e}",
             color=colors[j], alpha=0.8)
        

plt.grid(True, linestyle='--', alpha=0.3)
plt.xlabel("[C]_0")
plt.ylabel("[XY]_{eq}")

plt.legend()
    
plt.show()

