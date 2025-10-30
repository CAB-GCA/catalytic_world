import numpy as np
import math
from random import random
from scipy.integrate import solve_ivp


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
    reactions = []
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
    species = reactions[:, :-1]  # The last element is the ktype
    species = species[species != '']
    species = np.array(sorted(np.unique(species), key=len))
    return species


def calculate_m(reactions):
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
    m = calculate_m(reactions)
    n_species = len(species)
    c = np.zeros(shape=(m, n_species))

    for i in range(m):
        reaction = reactions[i, :]
        for j in range(len(reaction)):
            if reaction[j] in species:  # If the element is different from '' it will be in species
                # Get the index of the element in the list of species (to keep it coherent)
                index = list(species).index(reaction[j])
                if j in range(0, 2):  # If it is a reactant
                    c[i, index] -= 1
                elif j in range(2, 4):  # Or a product
                    c[i, index] += 1
    return c


def reactants(c):
    c_reactants = abs(np.where(c < 0, c, 0))
    return c_reactants


def chemistry(method, iterations, reactions, food_molecules, initial_food, k, V):
    species = obtain_species(reactions)
    c = c_matrix(reactions, species)
    k_types = reactions[:, -1]
    m = calculate_m(reactions)

    if m < len(k):
        raise Exception(
            'No se han definido suficientes constantes de reacción')

    if method == 'Gillespie':
        abundances = np.zeros((1, np.shape(species)[0]))
        abundances[0, :food_molecules] = initial_food
        abundances, times, V = gillespie(abundances, m, species,
                                      k_types, k, c, V, iterations)

    elif method == 'Deterministic':
        abundances = np.zeros((1, np.shape(species)[0]))
        abundances[0, :food_molecules] = initial_food
        times, abundances, V = integrate_ODEs(reactions, k, V, abundances[0, :],
                                           iterations, species)

    return abundances, times, V


def gillespie(abundances, m, species, k_types, k, c, V, iterations):
    # Get h (one per reaction)
    c_reactants = reactants(c)
    h = np.zeros(m)
    a = np.zeros(m)
    times = np.array([0.0], dtype=float)
    V = np.array([np.ravel(V)[-1]], dtype=float)
    abundance_v_relation = np.sum(abundances) / V

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
        
        elif k_types[i] == '4':
            #h_m = 1
            h[i] = V
            
        # Get the a_m
        a[i] = h[i]*k[i]

        return a

    def update_a(c):
        update = []
        # posiciones en las que c no sea 0, se cambian a 1
        c = np.where(c != 0, 1, 0)
        # esto evita que terminos queden como 0 por sumas (por ejemplo 2-2 = 0, aunque
        # sería posible obtener este resultado por casualidad)

        for m in c:  # por cada fila, cada reaccion
            # tengo que obtener las reacciones que contienen especies que se estan transformando:
            # non.zero devuelve qué posiciones son diferentes a 0
            update.append(np.nonzero(c @ m))

        return update
    
    def update_v(abundance, abundance_v_relation):
        total_abundance = np.sum(abundance)

        return total_abundance/abundance_v_relation
        

    # lista con tantos elementos como reacciones haya
    reactions_to_update = update_a(c)
    # si sale la reacción 0, el elemento 0 es un array con las reacciones que se tienen
    # que actualizar en ese caso, [0,2,3] indicará que se tienen que actualizar las reacciones
    # 0, 2 y 3, dejando sin actualizar la 1 porque no ha habido cambios en sus reactivos

    for n in range(iterations):
        abundance = abundances[n]

        if n != 0:
            for i in np.nditer(reactions_to_update[mu]):
                a = calculate_a(a, i, k_types, abundance, c_reactants, h, k, V[-1])

        else:
            for i in range(m):
                a = calculate_a(a, i, k_types, abundance, c_reactants, h, k, V[-1])

        # Get the a_0
        a0 = np.sum(a)
        if a0 == 0:
            print("La probabilidad total es 0 !!")
            return abundances, times, V

        # Get two random numbers, r1 and r2
        r1 = random()
        r2 = random()

        # Get mu and tau
        tau = (1/a0) * math.log(1/r1)

        sum_a = np.cumsum(a)
        for mu in range(len(a)):
            if sum_a[mu] >= r2*a0:
                break

        abundances = np.vstack((abundances, abundances[n] + c[mu]))
        times = np.append(times, times[-1] + tau)
        V = np.append(V, update_v(abundances[-1, :], abundance_v_relation))


    return abundances, times, V


def integrate_ODEs(reactions, k, V, initial_abundance, iterations, species):
    n_species = len(species)
    c = c_matrix(reactions, species)
    c_reactants = reactants(c)
    k_types = reactions[:, -1]
    m = calculate_m(reactions)
    t_end = iterations

    def calculate_dxdt(dxdt, i, k_types, abundance, c_reactants, h, k, V):
        if k_types[i] == '1':
            # h_m = x1
            x = abundance[c_reactants[i] == 1]
            h[i] = x/V

        elif k_types[i] == '2':
            # h_m = x^2
            x = abundance[c_reactants[i] == 2]
            h[i] = (x/V)**2

        elif k_types[i] == '3':
            # h_m = x1*x2
            x = abundance[c_reactants[i] == 1]
            h[i] = np.prod(x)*((1/V)**2)

        # Get the a_m
        dxdt[i] = h[i]*k[i]

        return dxdt

    def rhs(t, X):
        h = np.zeros(m)
        dxdt = np.zeros(m)
        for i in range(m):
            dxdt = calculate_dxdt(dxdt, i, k_types, X, c_reactants, h, k, V)
        return c.T @ dxdt

    sol = solve_ivp(rhs, [0, t_end], initial_abundance, method='BDF')
    times = sol.t
    abundances = sol.y.T


    return times, abundances, V
