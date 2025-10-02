import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from fun_gilles import *

# Initialization:
file = "reactions_XY.txt" # M reactions
n_iterations = 1000
method = "Gillespie"
# Reaction constants:
k = [0.1,0.01]
# Volume:
V = 100
initial_food = 50 # Initial population number
food_molecules = 2

reactions = read_file(file)
species = obtain_species(reactions)
abundances = np.zeros((n_iterations,np.shape(species)[0]))
abundances[0,:food_molecules] = initial_food
c = c_matrix(reactions, species)

print(gillespie(abundances, reactions, species, reactions[:,-1], k, n=0, t=0, c=c))