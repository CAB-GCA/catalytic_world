import numpy as np
import matplotlib.pyplot as plt
from fun_gilles import *

# Initialization:
file = "reactions_XYC_food_XY.txt" # M reactions
method = "Gillespie" # Gillespie or Deterministic

# Reaction constants:
k = [1,0,1,0]+ [1e-2] # len(k)= # de reacciones
# Volume:
V = 1000

# condiciones iniciales
initial_food = [100,0,0,0,0] # initial molecules number

# obtener reacciones y especies:
reactions = read_file(file)
species = obtain_species(reactions)
n_iterations= 10000 # In the deterministic mode n_iterations refers to the t_end


abundances, times, V = chemistry(method, n_iterations, reactions,
                                initial_food, k, V)

print(f"Parameters used for simulation:\n\
Initial concentrations:\nX_0={initial_food[1]}\nY_0={initial_food[2]}\n\
C_0={initial_food[0]}\n\
k_a = {k[0]}; k_a_r = {k[1]}\n\
k_b= {k[2]}, k_b_r={k[3]}\n\
# iterations = {n_iterations}")