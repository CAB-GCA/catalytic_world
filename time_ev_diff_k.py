import numpy as np
import matplotlib.pyplot as plt; plt.rcParams['text.usetex'] = True
from fun_gilles import *

# Initialization:
file = "reactions_XYC.txt" # M reactions
method = "Gillespie" # Gillespie or Deterministic

# Reaction constants:
k = [100,1e4,100,100] # len(k)= # de reacciones
# Volume:
V = 1

# condiciones iniciales
initial_food = [6000,10000,10000] # initial molecules number
food_molecules = 3

# obtener reacciones y especies:
reactions = read_file(file)
species = obtain_species(reactions)
n_iterations= 35600


abundances, times = chemistry(method, n_iterations, reactions, food_molecules, 
                              initial_food, k, V)

# Representation
plt.figure()
colors = plt.cm.coolwarm(np.linspace(0, 1, len(species)))
plt.grid()

for i in range(len(species)):
    plt.plot(times[:], abundances[:, i], label=species[i] + " low $k_{ar}$", 
              color=colors[i], alpha=0.9)
    
k = [100,1e7,100,100] # len(k)= # de reacciones
n_iterations= 90000

abundances, times = chemistry(method, n_iterations, reactions, food_molecules, 
                              initial_food, k, V)

for i in range(len(species)):
    plt.plot(times[:], abundances[:, i], label=species[i] + " high $k_{ar}$", 
              color=colors[i], alpha=0.9, linestyle= (0, (1, 10)))

plt.grid(True, linestyle='--', alpha=0.3)
plt.xlabel("Time")
plt.ylabel("Abundances")
plt.legend(loc="upper center", ncols=2, fontsize='small')
plt.title(f"{method} time evolution: {np.shape(times)[0]-1} iterations, k_a= {k[0]}")
    


plt.show()