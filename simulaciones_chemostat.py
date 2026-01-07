import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')  # allow importing from parent directory
from fun_gilles import *
import pickle
import pandas as pd


results = []
output_file = "gillespie_chemostat.pkl"
try: 
    with open(output_file, "ab") as file:
        
        abundances, times, volumes = chemistry(method= "Gillespie",
                                    iterations= 5e6,
                                    file="examples/reactions_autocat.txt",
                                    initial_food=[0]*4 + [6000] + [0]*3,
                                    k= [1]*8+[0.2]*4,
                                    V=1000,
                                    threshold= 0.)
        
        results.append((abundances,times,volumes))
        
        pickle.dump(results, file)

        file.flush()

        abundances, times, volumes = chemistry(method= "Gillespie",
                                            iterations= 5e6,
                                            file="examples/reactions_autocat.txt",
                                            initial_food=[0]*4 + [1000] + [0]*3,
                                            k= [1]*8+[0.2]*4,
                                            V=1000,
                                            threshold= 0.)
        
        results.append((abundances,times,volumes))
        
        pickle.dump(results, file)

        file.flush()
except Exception as e:
    print(f"An error has occured: {e}")

finally: 
    print("Simulation run completed")  
        

