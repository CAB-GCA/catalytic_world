from fun_gilles import *
import matplotlib.pyplot as plt
import numpy as np

file = "reactions_XYC_food_XY.txt" # M reactions
n_iterations = 200000
method = "Gillespie" # Gillespie or Deterministic
# Reaction constants:
k = [1]*8 # len(k)= # de reacciones
# Volume:
V = 2000
initial_food = [5000, 500, 5000, 0, 0] # initial molecules number

reactions = read_file(file)
species = obtain_species(reactions)

c = c_matrix(reactions, species)


# def update_a(c):
#     update = []
#     c = np.where(c != 0, 1, 0) # posiciones en las que c no sea 0, se cambian a 1
#     # esto evita que terminos queden como 0 por sumas (por ejemplo 2-2 = 0, aunque
#     # sería posible obtener este resultado por casualidad)
    
#     for m in c: # por cada fila, cada reaccion
#     # tengo que obtener las reacciones que contienen especies que se estan transformando:
#         update.append(np.nonzero(c @ m)) # non.zero devuelve tantos arrays como 
#         # dimensiones tenga el resultado, indicando qué posiciones son diferentes a 0
    
#     return update

# a = update_a(c)
# print(a)
# for i in np.nditer(a[3]):
#     print(i)

# def get_non_volume_species_indices(k_types, c):
#         """
#         Identifies the indices of species that are products in a type '4' reaction,
#         as these should be EXCLUDED from volume calculation.

#         Input
#         ----------
#         reactions : numpy array containing the reactions information.
#         species : the different species involved.
#         c : The stoichiometry matrix (reactions in rows, species in columns).

#         Returns
#         -------
#         non_volume_species_indices : numpy array of indices (columns of c) 
#                                     corresponding to the species whose abundance 
#                                     should NOT be used for volume calculation.
#         """
#         type_4_reaction_indices = np.where(k_types == '4')[0]
        
#         # Species indices that are products (c[i, j] > 0) in any type 4 reaction
#         non_volume_species_indices = []
#         for reaction_idx in type_4_reaction_indices:
#             # Find indices j where c[reaction_idx, j] > 0
#             product_indices = np.where(c[reaction_idx, :] > 0)[0]
#             non_volume_species_indices.append(product_indices)

#         # Return unique indices
#         return np.unique(non_volume_species_indices)
    
# print(get_non_volume_species_indices(reactions[:,-1], c))
# print(np.arange(np.shape(c)[1]))

abundances, times, V = chemistry(method= 'Gillespie', iterations= n_iterations,
                                 reactions= reactions, initial_food= initial_food,
                                 k = k, V = V)

colors = plt.cm.Spectral(np.linspace(0, 1, len(species)))

# --- Create the figure and a 1-row, 3-column subplot grid ---
# Using a wide figsize for a horizontal layout.
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 5), sharey=False) 
# Note: sharex=True is removed because we need the x-axis ticks on all three plots, 
# although they are the same (time).

# --- Subplot 1 (Left): Concentration ---
ax1 = axes[0]
for i in range(len(species)):
    # Calculate Concentration = Abundance / Volume
    ax1.plot(times, abundances[:, i] / V, label=species[i], color=colors[i], alpha=0.9)
ax1.grid(True, linestyle='--', alpha=0.3)
ax1.set_xlabel("Time")
ax1.set_ylabel("Concentration")
ax1.set_title("Concentration Evolution")

# --- Subplot 2 (Middle): Volume ---
ax2 = axes[1]
ax2.plot(times, V, color='gray')
ax2.grid(True, linestyle='--', alpha=0.3)
ax2.set_xlabel("Time")
ax2.set_ylabel("Volume")
ax2.set_title("Volume Evolution")


# --- Subplot 3 (Right): Abundances and Legend ---
ax3 = axes[2]
for i in range(len(species)):
    ax3.plot(times, abundances[:, i], label=species[i], color=colors[i], alpha=0.9)
ax3.grid(True, linestyle='--', alpha=0.3)
ax3.set_xlabel("Time")
ax3.set_ylabel("Abundances")
ax3.set_title("Abundance Evolution")
# Place the legend here, slightly outside the plot area for all species
ax3.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), title="Species")


# --- Final Touches ---
# Add a figure-wide title
fig.suptitle(f"Simulation Results (k= {k})", fontsize=16)

# Adjust layout to prevent labels/titles from overlapping
plt.tight_layout(rect=[0, 0, 1, 1]) # Adjust tight_layout for a single row

plt.show()