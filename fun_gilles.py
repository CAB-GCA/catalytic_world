import numpy as np
import math

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
        - 1: bimolecular (equal a+a -> )
        - 2: monomolecular (a -> )
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

def c_matrix(reactions):
    '''
    
    Parameters
    ----------
    reactions : numpy array containing the reactions in rows. The first and second
    column will be the reactants and the third and fourth the products. The last
    row indicates the type of reaction (monomolecular, bimolecular...)

    Returns
    -------
    c : a matrix with reactions in rows and a specie in each column. The values
    indicate the steichiometry of the specie in that reaction.

    '''
    n_reactions = calculate_n_reactions(reactions)
    species = obtain_species(reactions)
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

        
def obtain_abundances(abundances, c_matrix, reaction):
    
    

