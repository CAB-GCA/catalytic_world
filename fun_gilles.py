import numpy as np
import math
from random import random, seed
from scipy.integrate import solve_ivp
from numpy.linalg import matrix_rank

seed(1)

def read_file(file_name: str):
    '''

    Input
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

    Input
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

    Input
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

    Input
    ----------
    reactions : numpy array containing the reactions in rows. The first and second
    column will be the reactants and the third and fourth the products. The last
    row indicates the type of reaction (monomolecular, bimolecular...)
    species : TYPE
        DESCRIPTION.

    Returns
    -------
    c :  a matrix with reactions in rows and a specie in each column. The values
    indicate the stoichiometry of the specie in that reaction.

    '''
    m = calculate_m(reactions) # number of reactions
    n_species = len(species) # number of species
    c = np.zeros(shape=(m, n_species)) # c matrix: a row for each reaction and a column for each species
    
    for i in range(m):
        reaction = reactions[i, :] 
        for j in range(len(reaction)):
            # each j will be a species or empty ('')
            if reaction[j] in species:  # If the element is different from '' it will be in species 
                                        # this eliminates empty species
                # Get the index of the element in the list of species (to keep it coherent)
                index = list(species).index(reaction[j])
                if j in range(0, 2):  # If it is a reactant
                    c[i, index] -= 1
                elif j in range(2, 4):  # Or a product
                    c[i, index] += 1
    return c


def reactants(c):
    """
    Input
    ----------
    c: c matrix
    
    Returns
    -------
    c_reactants: a C matrix with only the reactant species for each
    reaction
    """
    c_reactants = abs(np.where(c < 0, c, 0))
    return c_reactants


def check_thermodynamics(reactions):
    c_matrix_wo_food = c_matrix(reactions[reactions[:,-1]!='4'], obtain_species(reactions))
    # max rank = # of species (which is the number of rows of the c_matrix)
    # since most simulations are run with reversible reactions (if not all)
    # the max rank = # of species / 2
    if matrix_rank(c_matrix_wo_food) >= c_matrix_wo_food.shape[0]/2:
        return True
    
    else:
        return False


def chemistry(method: str, iterations:int, reactions, 
              initial_food:list, k:list, V, **kwargs):
    """
    Performs the method indicated
    
    Input
    ----------
    method: "Gillespie" or "Protocell" for estochastic simulations or 
        "Deterministic"
    iterations: for stochastic simulations iterations will be the number of 
        performed iterations of the algorithm
        for deterministic simulations, iterations will be the final time
        of the simulation
    reactions: numpy array containing the reactions in rows. The first and second
        column will be the reactants and the third and fourth the products. The last
        row indicates the type of reaction (monomolecular, bimolecular...)
    initial_food: numpy array which contains the inital number of molecules
        this should be as long as the number of species of the simulation
    k: list containing all the catalytic constants for each reaction 
    V: initial volume
    
    Returns
    -------
    time: a 1-dimensional numpy array containing the times at which each iteration
        took place
    abundances: a n_species x iterations performed numpy array. Each row will be the
        abundance found of each species at that iteration
    V: a 1-dimensional numpy array containing volume of the system at each iteration
        The volume is defined so that the relation (total number of molecules)/(volume)
        remains constant
    """
    species = obtain_species(reactions)
    c = c_matrix(reactions, species)
    k_types = reactions[:, -1] # monomolecular, bimolecular...
    m = calculate_m(reactions) # number of reactions
    abundances = np.zeros((1, np.shape(species)[0])) # initialization of abundances
    
    if len(initial_food) != len(species): # not enough initial abundances
        raise Exception(
            'No se han definido el mismo numero de condiciones\
                iniciales que de especies'
        )
    
    abundances[0, :] = initial_food

    if check_thermodynamics(reactions) == False:
        raise Exception('El sistema no es compatible termodinámicamente')

    if method == 'Gillespie': # estochastic algorithm
        # TODO cuando funcione protocell con el threshold bien ->
        # - quitar iterations de protocell
        # - dejar iterations solo en Gillespie y deterministic
        # iterations = kwargs.get('iterations', None) 
        
        # if iterations is None:
        #     raise Exception("Gillespie method requires the 'iterations' parameter.")
        
        abundances, times, V = gillespie(abundances, m, k_types, k, c, V, iterations)

    elif method == 'Deterministic': # deterministic algorithm
        
        # iterations = kwargs.get('iterations', None) 
        # if iterations is None:
        #     raise Exception("Deterministic method requires the 'iterations' parameter.")
        
        times, abundances, V = integrate_ODEs(reactions, k, V, abundances,
                                            iterations, species)

    elif method == 'Protocell': # estochastic algorithm
        threshold = kwargs.get('threshold', None) 
        
        abundances, times, V = gillespieProtocell(abundances, m, k_types, k, c, 
                                                V, iterations, threshold)
    
    return abundances, times, V

def update_a(c):
    """
    This allows the algorithm to only update the probabilities of the
    reactions that have been affected in the previous iteration
    
    Input
    ------
    c: C matrix
    
    Returns
    -------
    update: list with len() = # of reactions. When a reaction is performed in 
        the previous iteration, this contains the reactions that need to be 
        updated (because the abundance of its species has change)
    """
    update = []
    # posiciones en las que c no sea 0, se cambian a 1
    c = np.where(c != 0, 1, 0)
    # esto evita que terminos queden como 0 por sumas (por ejemplo 2-2 = 0)

    for m in c:  # por cada fila, cada reaccion
        # tengo que obtener las reacciones que contienen especies que se estan transformando:
        # non.zero devuelve qué posiciones son diferentes a 0
        update.append(np.nonzero(c @ m))

    return update

def calculate_a(a, i, k_types, abundance, c_reactants, h, k, V):
    """
    Calculates the probability of every possible reaction
    given the quantity of each species, the volume, the reaction
    type and the catalytic constant
    
    Input
    -------
    (a, k_types, abundance, c_reactants, h, k and V: 
        Previously explained)
    
    i: iteration of a calculation (this function is used in a loop)
        will be equivalent to the reaction number
    
    Returns
    -------
    a: probability of each reaction taking place
    """
    if k_types[i] == '1':
        # h_m = X1
        x = abundance[c_reactants[i] == 1]
        h[i] = x

    elif k_types[i] == '2':
        # h_m = (1/2)*X1*(X1-1)
        x = abundance[c_reactants[i] == 2]
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

def gillespie(abundances, m, k_types, k, c, V, iterations):
    """
    Performs the Gillespie algorithm
    
    Input
    ----------
    abundances: the initial abundances for each species
    m: number of reactions
    k_types: reaction type
        - 1: monomolecular (a -> )
        - 2: bimolecular (equal a+a -> )
        - 3: bimolecular (but different a+b -> )
        - 4: food generation ( -> food)
    k: catalytic constant for each reaction
    c: C matrix, a matrix with reactions in rows and a specie in each column. The values
        indicate the stoichiometry of the specie in that reaction.
    V: initial volume
    iterations: number of iterations of the algorithm
    
    Returns
    -------
    abundances: numpy array with dimensions = # of species x # of iterations
        the quantity of each species at each iteration
    times: numpy array with dimensions = # of iterations x 1
        the times at which a reaction takes place according to the algorithm
    V: a 1-dimensional numpy array containing volume of the system at each iteration
        The volume is defined so that the relation (total number of molecules)/(volume)
        remains constant
    """
    if m != len(k): # not enough catalytic constants
        raise Exception(
            'No se ha definido el mismo numero de constantes\n\
            de reaccion que de reacciones'
            )
    
    # --- INITIALIZATION ---
    
    # Get h (one per reaction)
    c_reactants = reactants(c) # c matrix with only reactivos
    # inizialization:
    h = np.zeros(m) # h = propensities (no tiene en cuenta la constante catalitica)
    a = np.zeros(m) # a = probabilites (propensities * constante catalítica)
    times = np.array([0.0], dtype=float)
    V = np.array([float(V)], dtype=float)
    
    reactions_to_update = update_a(c)
    # si sale la reacción 0, el elemento 0 es un array con las reacciones que se tienen
    # que actualizar en ese caso, [0,2,3] indicará que se tienen que actualizar las reacciones
    # 0, 2 y 3, dejando sin actualizar la 1 porque no ha habido cambios en sus reactivos
    
    # --- VOLUME CALCULATION ---
    
    non_volume_species_indices = get_non_volume_species_indices(k_types, c)
    # Calculate initial abundance_v_relation based on the contributing species
    all_species_indices = np.arange(np.shape(c)[1])
    # The species that contribute are the ones that are not on non_volume_species
    # np.setdiff1d returns the different values in two arrays --> 
    # in this case the values will be the indices from species that affect the protocell volume
    volume_species_indices = np.setdiff1d(all_species_indices, non_volume_species_indices)
    
    if len(volume_species_indices) == 0:
        raise Exception("No species are products of a non type '4' reaction, so volume cannot be calculated.")
    
    initial_total_abundance = np.sum(abundances[0, volume_species_indices])
    abundance_v_relation = initial_total_abundance / V[0]

    # --- GILLESPIE ALGORITM ---

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

        abundances = np.vstack((abundances, abundance + c[mu]))
        times = np.append(times, times[-1] + tau)
        V = np.append(V, update_v_protocell(abundance, abundance_v_relation, volume_species_indices))


    return abundances, times, V

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

def integrate_ODEs(reactions, k, V, initial_abundance, iterations, species):
    n_species = len(species)
    c = c_matrix(reactions, species)
    c_reactants = reactants(c)
    k_types = reactions[:, -1]
    m = calculate_m(reactions)
    t_end = iterations
    
    if m != len(k): # not enough catalytic constants
        raise Exception(
            'No se ha definido el mismo numero de constantes\n\
            de reaccion que de reacciones'
        )
    
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

def calculate_a_without4(a, i, k_types, abundance, c_reactants, h, k, V):
    """
    Calculates the probability of every possible reaction
    given the quantity of each species, the volume, the reaction
    type and the catalytic constant
    
    Input
    -------
    (a, k_types, abundance, c_reactants, h, k and V: 
        Previously explained)
    
    i: iteration of a calculation (this function is used in a loop)
        will be equivalent to the reaction number
    
    Returns
    -------
    a: probability of each reaction taking place
    """
    if k_types[i] == '1':
        # h_m = X1
        x = abundance[c_reactants[i] == 1]
        h[i] = x

    elif k_types[i] == '2':
        # h_m = (1/2)*X1*(X1-1)
        x = abundance[c_reactants[i] == 2]
        h[i] = (1/2)*x*(x-1)*(1/V)

    elif k_types[i] == '3':
        # h_m = X1*X2
        x = abundance[c_reactants[i] == 1]
        h[i] = np.prod(x)*(1/V)
    
    # Get the a_m
    a[i] = h[i]*k[i]

    return a

def get_non_volume_species_indices(k_types, c):
    """
    Identifies the indices of species that are products in a type '4' reaction,
    as these should be EXCLUDED from volume calculation.

    Input
    ----------
    reactions : numpy array containing the reactions information.
    species : the different species involved.
    c : The stoichiometry matrix (reactions in rows, species in columns).

    Returns
    -------
    non_volume_species_indices : numpy array of indices (columns of c) 
                                corresponding to the species whose abundance 
                                should NOT be used for volume calculation.
    """
    type_4_reaction_indices = np.where(k_types == '4')[0]
    
    # Species indices that are products (c[i, j] > 0) in any type 4 reaction
    non_volume_species_indices = []
    for reaction_idx in type_4_reaction_indices:
        # Find indices j where c[reaction_idx, j] > 0
        product_indices = np.where(c[reaction_idx, :] > 0)[0]
        non_volume_species_indices.append(product_indices)

    # Return unique indices
    return np.unique(non_volume_species_indices)

def update_v_protocell(abundance, abundance_v_relation, volume_species_indices):
    """
    Updates the volume of the system based on ALL species EXCEPT those specified.

    Input
    ------
    abundance: the last row of the abundances matrix (all species).
    abundance_v_relation: the constant relation (Total Abundance / Volume)
        based on the species that DO affect the volume.
    volume_species_indices: indices of species to include from 
        total abundance calculation.
        
    Returns
    --------
    total_abundance/abundance_v_relation: the volume for the last iteration
    """
    
    # Calculate total abundance using only the contributing species indices
    total_abundance = np.sum(abundance[volume_species_indices])

    if total_abundance <= 0:
        return 1e-6 

    return total_abundance/abundance_v_relation

def block_statistics(abundances):
    """
    Calculates the mean and standard deviation for non-overlapping blocks.
    W is the block size (e.g., 100).
    """
    block_std = np.std(abundances, axis=0)
    
    return block_std

def threshold_function(V: float) -> float:
    """
    Predicts the required SD threshold for a new volume V_new.
    Threshold = (A / V_new^m) * (1 + margin)
    """
    A_fit, m_fit = 6.595, -0.987
    margin = 0.5
    if V <= 0:
        raise ValueError("Volume must be positive.")
        
    predicted_sigma = A_fit * (V) ** m_fit
    threshold = predicted_sigma * (1 + margin)
    return threshold

def gillespieProtocell(
    abundances, 
    m, 
    k_types, 
    k, 
    c, 
    V, 
    iterations, 
    threshold
) -> tuple:
    times = np.array([0.0], dtype=float)
    V = np.array([float(V)], dtype=float)
    
    if threshold == None: # if user does not define a specific threshold
        threshold = threshold_function(V)
    
    non_volume_species_indices = get_non_volume_species_indices(k_types, c)
    # Calculate initial abundance_v_relation based on the contributing species
    all_species_indices = np.arange(np.shape(c)[1])
    # The species that contribute are the ones that are not on non_volume_species
    # np.setdiff1d returns the different values in two arrays --> 
    # in this case the values will be the indices from species that affect the protocell volume
    volume_species_indices = np.setdiff1d(all_species_indices, non_volume_species_indices)
    
    if len(volume_species_indices) == 0:
        raise Exception("No species are products of a non type '4' reaction, so volume cannot be calculated.")
    
    initial_total_abundance = np.sum(abundances[0, volume_species_indices])
    abundance_v_relation = initial_total_abundance / V[0]
    
    # Calculate initial concentration for food, as this will be constant throughout the simulation
    if non_volume_species_indices.any != None:
        initial_non_volume_conc = abundances[0, non_volume_species_indices] / V[0]      
    
    # now we filter the reactions -> keeping only those that are not type 4
    no_type_4_reaction_indices = np.where(k_types != '4')[0]
    reactions_to_keep = np.zeros(m, dtype= bool)
    reactions_to_keep[no_type_4_reaction_indices]= True
    
    # variables to update:
    k_types = k_types[reactions_to_keep]
    c = c[reactions_to_keep]
    m = np.shape(c)[0] # number of reactions must be updated
    c_reactants = reactants(c) # c matrix with only reactivos
    
    # inizialization of gillespie algorithm
    h = np.zeros(m) # h = propensities (no tiene en cuenta la constante catalitica)
    a = np.zeros(m) # a = probabilites (propensities * constante catalítica)
    reactions_to_update = update_a(c)

    if len(no_type_4_reaction_indices) != len(k):
        raise Exception(
            'El número de constantes de reacción es diferente al número de reacciones')
    
    n = 0
    counter = 0
    while n < iterations:
        abundance = abundances[n]

        if n != 0:
            for i in np.nditer(reactions_to_update[mu]):
                a = calculate_a_without4(a, i, k_types, abundance, c_reactants, h, k, V[-1])

        else:
            for i in range(m):
                a = calculate_a_without4(a, i, k_types, abundance, c_reactants, h, k, V[-1])

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
        
        new_abundances = abundances[n] + c[mu]
        times = np.append(times, times[-1] + tau)
        # actualizar el volumen en función de la reacción que haya tocado
        new_V = update_v_protocell(new_abundances, abundance_v_relation, volume_species_indices)
        V = np.append(V, new_V)
        
        if non_volume_species_indices.any != None:
            # now there's more/less food that can be inside that volume
            new_abundances[non_volume_species_indices] = np.round(initial_non_volume_conc * new_V)
            
        abundances = np.vstack((abundances, new_abundances))

        # stop criterion
        if n%500 == 0 and n > 2000:
            last_500_concentrations = (abundances[-500:,volume_species_indices].T/V[-500:]).T
            std = block_statistics(last_500_concentrations)
            
            if max(std) < threshold:
                counter += 1
                if counter == V[0]/10:
                    return abundances, times, V
            else: 
                counter = 0
        
        n += 1

    print("Criterion for stop was # of iterations")
    return abundances, times, V