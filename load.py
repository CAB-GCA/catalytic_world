from fun_gilles import *
import matplotlib.pyplot as plt
import numpy as np
import pickle

def plot(abundances, times, V):
    colors = plt.cm.rainbow(np.linspace(0, 1, len(species)))

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
    fig.suptitle(f"Simulation Results", fontsize=16)

    # Adjust layout to prevent labels/titles from overlapping
    plt.tight_layout(rect=[0, 0, 1, 1]) # Adjust tight_layout for a single row

    plt.show()

file = 'reactions_autocat.txt'
species = obtain_species(read_file(file))

data_file = open('example','rb')
data = pickle.load(data_file)

V = data[:,-1]
times = data[:,-2]
abundances = data[:,:-2]
print(data.shape)


data_file.close()

plot(abundances, times, V)