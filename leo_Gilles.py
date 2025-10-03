import numpy as np
import matplotlib.pyplot as plt
from fun_gilles import *

# Initialization:
file = "reactions_XYC.txt" # M reactions
n_iterations = 10000
method = "Gillespie"
# Reaction constants:
k = [0.01,0.015]
# Volume:
V = 100
initial_food = 50 # Initial population number
food_molecules = 3

reactions = read_file(file)
species = obtain_species(reactions)
abundances = np.zeros((n_iterations,np.shape(species)[0]))
abundances[0,:food_molecules] = initial_food
c = c_matrix(reactions, species)
times = np.zeros(n_iterations)
t = 0
mu, a, h = 0, 0, 0
for n in range(0, n_iterations-1):
    abundances, n, t, mu,a, h = gillespie(abundances, reactions, species, reactions[:,-1], k, n, t, c, mu, a, h)
    times[n] = t
    
print(abundances, times)

# Representation
plt.figure()
colors = plt.cm.viridis(np.linspace(0, 1, len(species)))
plt.grid()

for i in range(len(species)):
    plt.plot(times, abundances[:, i], label=species[i], color=colors[i], alpha=0.8)

plt.legend()
    
plt.show()