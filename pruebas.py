# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 11:12:09 2025

para entender por qué

@author: cvzad
"""

from fun_gilles import *

file = "C:/Users/cvzad/Documents/all_obs/master/tfm/catalytic_world/reactions_0.txt" # M reactions
n_iterations = 1000
method = "Gillespie" # Gillespie or Deterministic
# Reaction constants:
k = [1,0.1] # len(k)= # de reacciones
# Volume:
V = 1
initial_food = [1000, 900] # initial molecules number
food_molecules = 2

reactions = read_file(file)
species = obtain_species(reactions)

c = c_matrix(reactions, species)


def update_a(c):
    update = []
    c = np.where(c != 0, 1, 0) # posiciones en las que c no sea 0, se cambian a 1
    # esto evita que terminos queden como 0 por sumas (por ejemplo 2-2 = 0, aunque
    # sería posible obtener este resultado por casualidad)
    
    for m in c: # por cada fila, cada reaccion
    # tengo que obtener las reacciones que contienen especies que se estan transformando:
        update.append(np.nonzero(c @ m)) # non.zero devuelve tantos arrays como 
        # dimensiones tenga el resultado, indicando qué posiciones son diferentes a 0
    
    return update

a = update_a(c)
print(a)
for i in np.nditer(a[3]):
    print(i)