import numpy as np
import matplotlib.pyplot as plt; plt.rcParams['text.usetex'] = True
from fun_gilles import *

# Initialization:
file = "reactions_XYC.txt" # M reactions
method = "Gillespie" # Gillespie or Deterministic

# Reaction constants:
k = [100]*4 # len(k)= # de reacciones
# Volume:
V = 1

# condiciones iniciales
initial_food = [0,1000,1000] # initial molecules number
food_molecules = 3

# obtener reacciones y especies:
reactions = read_file(file)
species = obtain_species(reactions)

# Different values for initial C concentration and k_a:
initial_c = np.round(np.linspace(1,2500,50))
k_a = np.logspace(-2, 5, 8)
n_iterations = 4000

mean_equilibrium = np.zeros((len(initial_c)))
colors = plt.cm.Spectral(np.linspace(0, 1, len(k_a)))

for j in range(len(k_a)):
    k[0] = k_a[j]
    
    for i in range(len(initial_c)):
        initial_food[0] = initial_c[i]
        abundances, times = chemistry(method, n_iterations, reactions, 
                                      food_molecules, initial_food, k, V)
        mean_equilibrium[i] = np.mean(abundances[-min(50, len(abundances)):, -1])
    
    plt.plot(initial_c, mean_equilibrium, label=f"$k_a$={k_a[j]:.1e}",
             color=colors[j], alpha=0.8)
        

plt.grid(True, linestyle='--', alpha=0.3)
plt.xlabel("$C_0$")
plt.ylabel("$XY_{eq}$")

plt.legend()
    
plt.show()

