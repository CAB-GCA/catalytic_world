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
        f = "examples/reactions_XYC_food.txt"
        abundances, times, V = chemistry(method= "Gillespie", iterations= int(1e7),
                                         file= f,
                                         initial_food= [500,0,0,0,0], 
                                         k= [1,0,1,1]+[0.1]*2 + [0.4]*2, 
                                         V= 500, 
                                         threshold= 0)
        
        results.append((abundances,times, V))
        
        pickle.dump(results, file)

        file.flush()
except Exception as e:
    print(f"An error has occured: {e}")

finally: 
    print("Simulation run completed")  
        



