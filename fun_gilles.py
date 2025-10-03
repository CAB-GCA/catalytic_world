import numpy as np
import math
from random import random

def read_file(file_name: str):
    '''
    
    Parameters
    ----------
    file_name : str
        a file containing the information of the reactions. This file must contain
        a header. 
        the reactions in rows. The first and second column will be the reactants 
        and the third and fourth the products. The last row indicates the type 
        of reaction 
        - 1: monomolecular (a -> )
        - 2: bimolecular (equal a+a -> )
        - 3: bimolecular (but different a+b -> )
        - 4: food generation ( -> food)
        

    Returns
    -------
    reactions : numpy array containing the reactions in rows. The first and second
    column will be the reactants and the third and fourth the products. The last
    row indicates the type of reaction (monomolecular, bimolecular...)

    '''
    reactions=[]
    with open(file_name, "r") as f:
        r = f.readlines()
        for reaction in range(len(r)):
            reactions.append(r[reaction].strip().split(","))
    reactions = np.array(reactions[1:])

    return reactions

def obtain_species(reactions):
    '''

    Parameters
    ----------
    reactions : numpy array containing the reactions in rows. The first and second
    column will be the reactants and the third and fourth the products. The last
    row indicates the type of reaction (monomolecular, bimolecular...)

    Returns
    -------
    species : the different species that are involved in those reactions

    '''
    species = reactions[:,:-1] # The last element is the ktype
    species = species[species != '']
    species = np.array(sorted(np.unique(species), key= len))
    return species

def length(x):
    return len(x)

def calculate_n_reactions(reactions):
    '''

    Parameters
    ----------
    reactions : numpy array containing the reactions in rows. The first and second
    column will be the reactants and the third and fourth the products. The last
    row indicates the type of reaction (monomolecular, bimolecular...)

    Returns
    -------
    The number of reactions

    '''
    return np.shape(reactions)[0]

def c_matrix(reactions, species):
    '''

    Parameters
    ----------
    reactions : numpy array containing the reactions in rows. The first and second
    column will be the reactants and the third and fourth the products. The last
    row indicates the type of reaction (monomolecular, bimolecular...)
    species : TYPE
        DESCRIPTION.

    Returns
    -------
    c :  a matrix with reactions in rows and a specie in each column. The values
    indicate the steichiometry of the specie in that reaction.

    '''
    n_reactions = calculate_n_reactions(reactions)
    n_species = len(species)
    c = np.zeros(shape= (n_reactions, n_species))

    for i in range(n_reactions):
        reaction = reactions[i,:]
        for j in range(len(reaction)):
            if reaction[j] in species: # If the element is different from '' it will be in species
                # Get the index of the element in the list of species (to keep it coherent)
                index = list(species).index(reaction[j]) 
                if j in range(0,2): # If it is a reactant
                    c[i, index] -= 1
                elif j in range(2,4): # Or a product
                    c[i, index] += 1
    return c

        
def reactants(c):
    c_reactants = abs(np.where(c<0, c, 0))
    return c_reactants

def calculate_a(a, i, k_types, abundance, c_reactants, h, k, V):
    if k_types[i] == '1':
        # h_m = X1
        x = abundance[c_reactants[i] == 1]
        h[i] = x
     
    elif k_types[i] == '2':
        # h_m = (1/2)*X1*(X1-1)
        x = abundance[c_reactants[i] == 1]
        h[i] = (1/2)*x*(x-1)*(1/V)
    
    elif k_types[i] == '3':
        # h_m = X1*X2
        x = abundance[c_reactants[i] == 1]
        h[i] = np.prod(x)*(1/V)
    
    # Get the a_m
    a[i] = h[i]*k[i]
    
    return a
    
def chemistry(method, iterations, reactions, food_molecules, initial_food, k, V):
    species = obtain_species(reactions)
    abundances = np.zeros((iterations+1,np.shape(species)[0]))
    abundances[0,:food_molecules] = initial_food
    c = c_matrix(reactions, species)
    times = np.zeros(iterations+1)
    k_types = reactions[:,-1]
    n_reactions = calculate_n_reactions(reactions)
    if method == 'Gillespie':
        t=0        
        for n in range(iterations):
            abundances, n, t = gillespie(abundances, n_reactions, species, 
                                                   k_types, k, n, t, c, V)
            times[n] = t
        
    return abundances, times

    
    
def gillespie(abundances, m, species, k_types, k, n_iteration, t, c_matrix, V):
    # Get h (one per reaction)
    c_reactants = reactants(c_matrix)
    abundance = abundances[n_iteration]
    h = np.zeros(m)
    a = np.zeros(m)
    
    for i in range(m):
        a = calculate_a(a, i, k_types, abundance, c_reactants, h, k, V)

    # Get the a_0
    a0 = np.sum(a)
    if a0 == 0:
        raise Exception("La probabilidad total es 0 !!")
    
    # Get two random numbers, r1 and r2
    r1 = random()
    r2 = random()
    
    # Get mu and tau
    tau = (1/a0) * math.log(1/r1)
    
    sum_a = np.cumsum(a)
    for mu in range(len(a)):
        if sum_a[mu] >= r2*a0:
            break
    
    abundances[n_iteration+1] = abundances[n_iteration] + c_matrix[mu]
    t += tau
    n_iteration += 1
    
    return abundances, n_iteration, t
    