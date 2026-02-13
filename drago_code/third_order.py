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
nombre_archivo = "third_order" + nombre_archivo + ".pkl"


third_order = "examples/third_order.txt"
intermediates = "examples/reactions_autocat.txt"

# simulation
method = "Protocell"
initial_species = [1000]*4 + [200]+[0]*3
initial_species_third = [1000]*4 + [200] + [0]
k_third_order = [1]*4
k_intermediates = [1]*8
V = 1
iterations = 1e5
sp_intermediates = obtain_species(read_file(intermediates))
sp_thirdorder = obtain_species(read_file(third_order))

try:
    with open(nombre_archivo, "ab") as file:
        
        for i in np.logspace(-5,3,9): 
            k_third_order[1]= i
            # k_third_order[3]= i
            a, t, v = chemistry(method,
                                iterations,
                                third_order,
                                initial_species_third,
                                k_third_order,
                                V,
                                threshold= 0)
            results = {i: (a,t,v)}
            pickle.dump(results, file)
            print(f"Simulation for k_2 and k_4 = {i} completed and saved.")
    
except Exception as e:
    print(f"An error has occured: {e}")
    
finally:
    print("Simulations completed")
    

