import numpy as np
import matplotlib.pyplot as plt; plt.rcParams['text.usetex'] = True
from fun_gilles import *

# Initialization:
file = "reactions_XYC.txt" # M reactions
method = "Deterministic" # Gillespie or Deterministic

food_molecules = 3
# obtener reacciones y especies:
reactions = read_file(file)
species = obtain_species(reactions)

k = [1,0,1,1] # len(k)= # de reacciones
# Volume:
initial_food = [0,1000,1000] # initial molecules number

initial_c = np.round(np.linspace(1,250,10))
k_var = [1e-2, 1]
n_iterations = 1000
method = "Gillespie"

equilibrium = np.zeros((len(initial_c)))
colors = plt.cm.Spectral(np.linspace(0, 1, len(k_var)))

for j in range(len(k_var)):
    k[3] = k_var[j] 
    
    for i in range(len(initial_c)):
        initial_food[0] = initial_c[i]
        abundances, times, V = chemistry(method, n_iterations, reactions, 
                                    food_molecules, initial_food, k, V0)
        equilibrium[i] = abundances[-1, -1]/V0
    
    plt.plot(initial_c/V0, equilibrium, label="$k_{br}$"+f"={k_var[j]:.1e}",
            color=colors[j], alpha=0.8)
        

plt.grid(True, linestyle='--', alpha=0.3)
plt.xlabel("$[C]_0$")
plt.ylabel("$[XY]_{eq}$")

plt.legend()
    
plt.show()


print(f"Parameters used for simulation:\n\
Initial concentrations:\nX_0={initial_food[1]}\nY_0={initial_food[2]}\n\
k_a = {k[0]}; k_a_r = {k_var}\n\
k_b= {k[2]}, k_b_r={k[3]}\n\
# iterations = {n_iterations}")
