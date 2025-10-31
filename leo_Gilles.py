import numpy as np
import matplotlib.pyplot as plt
from fun_gilles import *

# Initialization:
file = "reactions_XYC.txt" # M reactions
method = "Gillespie" # Gillespie or Deterministic

reactions = read_file(file)

k = [1,1,1,1] # len(k)= # de reacciones
# Volume:
V = 1000
initial_food = [0,1000,1000] + [0]*2 # initial molecules number
initial_c = np.round(np.linspace(1,250,20))
k_var = [1e-2,1]
n_iterations = 40000
method = "Gillespie"


abundances, times, V = chemistry(method, n_iterations, reactions,
                                initial_food, k, V)

print(f"Parameters used for simulation:\n\
Initial concentrations:\nX_0={initial_food[1]}\nY_0={initial_food[2]}\n\
C_0={initial_food[0]}\n\
k_a = {k[0]}; k_a_r = {k[1]}\n\
k_b= {k[2]}, k_b_r={k[3]}\n\
# iterations = {n_iterations}")