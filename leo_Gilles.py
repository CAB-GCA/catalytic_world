import numpy as np
import matplotlib.pyplot as plt
from fun_gilles import *

# Initialization:
file = "reactions_XYC.txt" # M reactions
method = "Deterministic" # Gillespie or Deterministic

# Reaction constants:
k = [1,1,1,1e4] # len(k)= # de reacciones
# Volume:
V = 1

# condiciones iniciales
initial_food = [1,10,10] # initial molecules number
food_molecules = 3

# obtener reacciones y especies:
reactions = read_file(file)
species = obtain_species(reactions)
n_iterations= 1000000 # In the deterministic mode n_iterations refers to the t_end


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
plt.title(f"{method} time evolution: {np.shape(times)[0]-1} iterations, k= {k}")
    
print(f"Parameters used for simulation:\n\
Initial concentrations:\nX_0={initial_food[1]}\nY_0={initial_food[2]}\n\
C_0={initial_food[0]}\n\
k_a = {k[0]}; k_a_r = {k[1]}\n\
k_b= {k[2]}, k_b_r={k[3]}\n\
# iterations = {n_iterations}")

plt.show()