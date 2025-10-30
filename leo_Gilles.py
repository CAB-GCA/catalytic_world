import numpy as np
import matplotlib.pyplot as plt
from fun_gilles import *

# Initialization:
file = "reactions_XYC_food.txt" # M reactions
method = "Gillespie" # Gillespie or Deterministic

# Reaction constants:
k = [0.1]*2+[1,0,1,1] # len(k)= # de reacciones
# Volume:
V = 1000

# condiciones iniciales
initial_food = [200,1000,1000,0,0] # initial molecules number

# obtener reacciones y especies:
reactions = read_file(file)
species = obtain_species(reactions)
n_iterations= 100000 # In the deterministic mode n_iterations refers to the t_end


abundances, times, V = chemistry(method, n_iterations, reactions,
                                initial_food, k, V)

# Representation
plt.figure()
colors = plt.cm.coolwarm(np.linspace(0, 1, len(species)))

for i in range(len(species)):
    plt.plot(times, abundances[:, i], label=species[i], color=colors[i], alpha=0.9)

plt.grid(True, linestyle='--', alpha=0.3)
plt.xlabel("Time")
plt.ylabel("Abundances")
plt.legend()
plt.title(f"{method} time evolution: {np.shape(times)[0]-1} iterations, k= {k}")
plt.show()

plt.figure()
plt.grid(True, linestyle='--', alpha=0.3)
plt.plot(times, V, color= colors[0])
plt.xlabel("Time")
plt.ylabel("Volume")
plt.show()

plt.figure()
plt.grid(True, linestyle='--', alpha=0.3)
for i in range(len(species)):
    plt.plot(times, abundances[:, i]/V, label=species[i], color=colors[i], alpha=0.9)

plt.legend()
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.show()


print(f"Parameters used for simulation:\n\
Initial concentrations:\nX_0={initial_food[1]}\nY_0={initial_food[2]}\n\
C_0={initial_food[0]}\n\
k_a = {k[0]}; k_a_r = {k[1]}\n\
k_b= {k[2]}, k_b_r={k[3]}\n\
# iterations = {n_iterations}")