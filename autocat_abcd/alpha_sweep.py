import pickle
import numpy as np
from numpy.polynomial.polynomial import polyfit


def load_two_parameter_sweep(filename):
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
    Applies your logic to get the slope (m) which is alpha.
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

def visualize_alpha_sweep(raw_data, k_values):
     
    n = len(k_values)
    alpha_grid = np.zeros((n, n))
    
    # 2. Process data
    for entry in raw_data:
        for idx, replicates in entry.items():
            row = idx // n
            col = idx % n
            
            if row >= n or col >= n:
                print(f"Warning: Index {idx} (row {row}, col {col}) is out of grid bounds {n}x{n}")
                continue
            
            replicate_alphas = []
            for rep in replicates:
                # rep = (abundances, times, volumes)
                times = np.array(rep[1])[int(len(np.array(rep[1]))/5):]
                volumes = np.array(rep[2])[int(len(np.array(rep[2]))/5):]
                
                # Calculate alpha for this specific run
                try:
                    a = extract_alpha_from_replicate(times, volumes)
                    replicate_alphas.append(a)
                except:
                    replicate_alphas.append(0) # Handle failed runs
            
            # Store the mean growth rate for this (k_in, k3) pair
            alpha_grid[row, col] = np.mean(replicate_alphas)

    return alpha_grid


k_change = np.logspace(-5,7,25)

for i in [0, 1]:
    data = load_two_parameter_sweep(f"autocat_abcd/barrido_k{i}_k{i+2}.pkl")
    
    # alpha_sweep = visualize_alpha_sweep(data, k_change)
    # with open(f"autocat_abcd/alpha_sweep_k{i}_k{i+2}.pkl", "ab") as f:

    #     try:
    #         pickle.dump(alpha_sweep, f)
    #     except Exception as e:
    #         print(f"An error has occured: {e}")

    #     finally: 
    #         print("Simulation run completed")
            
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
    data = None
        
    


