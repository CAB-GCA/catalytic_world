import numpy as np
import matplotlib.pyplot as plt; plt.rcParams['text.usetex'] = True
from fun_gilles import *

# Initialization:
file = "reactions_XYC.txt" # M reactions
method = "Deterministic" # Gillespie or Deterministic

# Reaction constants:
k = [1,0,1,1] # len(k)= # de reacciones
# Volume:
V = 1

# condiciones iniciales
initial_food = [0,10000,10000] # initial molecules number
food_molecules = 3

# obtener reacciones y especies:
reactions = read_file(file)
species = obtain_species(reactions)

# Different values for initial C concentration and k_var:
initial_c = np.round(np.linspace(1,1000,6))
k_var = np.logspace(-4, 5, 30)
n_iterations = 100000000

equilibrium = np.zeros((len(k_var)))
colors = plt.cm.Spectral(np.linspace(0, 1, len(initial_c)))

    
for i in range(len(initial_c)):
    initial_food[0] = initial_c[i]
    
    for j in range(len(k_var)):
        k[3] = k_var[j]
        abundances, times = chemistry(method, n_iterations, reactions, 
                                      food_molecules, initial_food, k, V)
        equilibrium[j] = abundances[-1, -1]
    
    plt.plot(k_var, equilibrium, label="$C_0$"+f"={initial_c[i]}",
             color=colors[i], alpha=0.8)
    
plt.xscale("log")
plt.grid(True, linestyle='--', alpha=0.3)
plt.xlabel("$k_{br}$")
plt.ylabel("$XY_{eq}$")

plt.legend()
    
plt.show()
