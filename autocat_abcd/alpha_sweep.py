import pickle
import numpy as np
from numpy.polynomial.polynomial import polyfit

"""
reads the barrido files (created with barrido_i_i+2.py files)
extracts the alpha value
saves it into a file for each pair of ki - ki+1 parameters

this code also saves a copy of the barrido_i_i+1.pkl file with only
one repetition to avoid memory problems
"""


def load_two_parameter_sweep(filename):
    # read the two parameter sweep file
    raw_data = []
    with open(filename, "rb") as f:
        while True:
            try:
                raw_data.append(pickle.load(f))
            except EOFError:
                break
    return raw_data

def extract_alpha_from_replicate(times, volumes):
    """
    Applies the logic to get the slope (m) which is the alpha of the exponential
    """
    # Ensure inputs are numpy arrays and handle log(0) safety
    times = np.array(times)
    volumes = np.array(volumes)
    
    # Linearize: log(V) = alpha * t + log(A)
    log_v = np.log(volumes)
    
    # polyfit returns [intercept, slope] for degree 1
    coeffs = polyfit(times, log_v, 1)
    alpha = coeffs[1] 
    return alpha

def create_alpha_grid(raw_data, k_values):
     
    n = len(k_values)
    alpha_grid = np.zeros((n, n))
    
    # Process data
    for entry in raw_data:
        
        for idx, replicates in entry.items():
            row = idx // n # the row indicates the value of k_i (1st row = 1st value of k_i)
            col = idx % n # the col indicates the value of k_i+2
            
            if row >= n or col >= n:
                print(f"Warning: Index {idx} (row {row}, col {col}) is out of grid bounds {n}x{n}")
                continue
            
            replicate_alphas = []
            for rep in replicates:
                # rep = (abundances, times, volumes)
                # we take the last 4/5 of the values to fit the alpha
                times = np.array(rep[1])[int(len(np.array(rep[1]))/5):] 
                volumes = np.array(rep[2])[int(len(np.array(rep[2]))/5):]
                
                # calculate alpha for this specific run
                try:
                    a = extract_alpha_from_replicate(times, volumes)
                    replicate_alphas.append(a)
                except:
                    replicate_alphas.append(0) # handle failed runs
            
            # store the mean growth rate for this (k_in, k3) pair
            alpha_grid[row, col] = np.mean(replicate_alphas)

    return alpha_grid


k_change = np.logspace(-5,7,25) # values tested for each k

for i in [0, 1]: # the two barridos
    data = load_two_parameter_sweep(f"autocat_abcd/barrido_k{i}_k{i+2}.pkl")
    # create the alpha grid for the file and save it
    alpha_sweep = create_alpha_grid(data, k_change)
    with open(f"autocat_abcd/alpha_sweep_k{i}_k{i+2}.pkl", "ab") as f:

        try:
            pickle.dump(alpha_sweep, f)
        except Exception as e:
            print(f"An error has occured: {e}")

        finally: 
            print("Simulation run completed")
    # create the short version of the barrido (by taking just one rep)        
    with open(f"autocat_abcd/short_barrido_k{i}_k{i+2}.pkl", "ab") as f:
        for element in data:
            for cond in element.keys():
                rep = element[cond][0]
                try:
                    pickle.dump({cond:rep}, f)
                except Exception as e:
                    print(f"An error has occured: {e}")

                finally: 
                    print("Simulation run completed")
    data = None # restart the data variable to prevent memory errors
        
    


