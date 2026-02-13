from datetime import datetime
import numpy as np
import sys
import os
import pickle
# Get the directory where your current script is located
current_dir = os.path.dirname(os.path.abspath(__file__))

# Move up one level to 'catalytic_world'
parent_dir = os.path.dirname(current_dir)

# Add it to the path
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from fun_gilles import *



nombre_archivo = "_".join(str(datetime.now()).split()).split(":") # ":" character is not accepted in a file name
nombre_archivo = "".join(nombre_archivo[:-1])
nombre_archivo = "simulation" + nombre_archivo + ".pkl"


intermediates = "./examples/reactions_autocat.txt"

# simulation
method = "Protocell"
initial_species = [100]*4 + [50]*2 +[0]*2
k = [1]*8
k[0] = 1e-4
k[2] = 1e6
V = 100
iterations = 1e6
sp_intermediates = obtain_species(read_file(intermediates))

try:
    with open(nombre_archivo, "ab") as file:

        a, t, v = chemistry(method,
                            iterations,
                            intermediates,
                            initial_species,
                            k,
                            V,
                            threshold= 0)
        results = (a,t,v)
        pickle.dump(results, file)
        print(f"Simulation completed and saved.")
    
except Exception as e:
    print(f"An error has occured: {e}")
    
finally:
    print("Simulations completed")
    

