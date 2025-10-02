import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from fun_gilles import *

# Initialization:
file = "reactions_XYC.txt" # M reactions
n_iterations = 10
method = "Gillespie"
# Reaction constants:
k = [0.01,0.1]
# Volume:
V = 100
initial_food = 50 # Initial population number
food_molecules = 3

reactions = read_file(file)
species = obtain_species(reactions)
abundances = np.zeros((n_iterations,np.shape(species)[0]))
abundances[0,:food_molecules] = initial_food
c = c_matrix(reactions, species)

for n in range(n_iterations-1):
    abundances, n, t = gillespie(abundances, reactions, species, reactions[:,-1], k, n, t=0, c=c)

print(abundances)