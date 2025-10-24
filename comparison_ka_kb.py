import numpy as np
import matplotlib.pyplot as plt; plt.rcParams['text.usetex'] = True
from fun_gilles import *

# Initialization:
file = "reactions_XYC.txt" # M reactions
method = "Deterministic" # Gillespie or Deterministic

# Reaction constants:
k = [100,100,100,100] # len(k)= # de reacciones
# Volume:
V = 1

# condiciones iniciales
initial_food = [0,10000,10000] # initial molecules number
food_molecules = 3

# obtener reacciones y especies:
reactions = read_file(file)
species = obtain_species(reactions)

# Different values for initial C concentration and k_var:
initial_c = np.round(np.linspace(1,25000,100))
k_arvar = np.logspace(4, 8, 3)
k_brvar = np.logspace(-2, 2, 3)
n_iterations = 100000

equilibrium = np.zeros((len(initial_c)))
colors = plt.cm.Spectral(np.linspace(0, 1, len(k_brvar)+3))
linestyles = ['-', '--', ':', '-.', (0, (3, 5, 1, 5))] 

for j in range(len(k_brvar)):
    k[3] = k_brvar[j]
    for h in range(len(k_arvar)):
        k[1] = k_arvar[h]
        for i in range(len(initial_c)):
            initial_food[0] = initial_c[i]
            abundances, times = chemistry(method, n_iterations, reactions, 
                                          food_molecules, initial_food, k, V)
            equilibrium[i] = abundances[-1, -1]
        
        plt.plot(initial_c, equilibrium, 
                 label="$k_{ar}$"+f"={k_arvar[h]:.1e}\n"+"$k_{br}$"+f"={k_brvar[j]:.1e}",
                 color=colors[j], linestyle= linestyles[h % len(linestyles)],
                 alpha=0.8)
        

plt.grid(True, linestyle='--', alpha=0.3)
plt.xlabel("$C_0$")
plt.ylabel("$XY_{eq}$")

plt.legend(fontsize='small', loc= "upper right")
    
plt.show()


print(f"Parameters used for simulation:\n\
Initial concentrations:\nX_0={initial_food[1]}\nY_0={initial_food[2]}\n\
k_a = {k[0]}; k_a_r = {k_arvar}\n\
k_b= {k[2]}, k_b_r={k_brvar}\n\
# iterations = {n_iterations}")
