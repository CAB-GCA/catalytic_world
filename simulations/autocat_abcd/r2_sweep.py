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

def generate_r2_grid(raw_data, k_values):
    """
    Processes raw stochastic data to calculate a mean R^2 grid.
    """
    n = len(k_values)
    r2_grid = np.zeros((n, n))
    
    for entry in raw_data:
        for idx, replicates in entry.items():
            row = idx // n
            col = idx % n
            
            if row >= n or col >= n:
                continue
                
            replicate_r2s = []
            for rep in replicates:
                # Use the same 20% burn-in cut as your alpha function
                times = np.array(rep[1])[int(len(np.array(rep[1]))/5):]
                volumes = np.array(rep[2])[int(len(np.array(rep[2]))/5):]
                
                try:
                    # Linearize: log(V) = alpha*t + log(A)
                    log_v = np.log(volumes)
                    coeffs = np.polyfit(times, log_v, 1)
                    
                    # Manual R^2 calculation
                    log_v_pred = coeffs[1] + coeffs[0] * times
                    ss_res = np.sum((log_v - log_v_pred)**2)
                    ss_tot = np.sum((log_v - np.mean(log_v))**2)
                    r2 = 1 - (ss_res / ss_tot)
                    replicate_r2s.append(r2)
                except:
                    replicate_r2s.append(0)
            
            # Store mean R^2, clipped at 0
            r2_grid[row, col] = max(0, np.mean(replicate_r2s))
            
    return r2_grid


k_change = np.logspace(-5,7,25)

for i in [0, 1]:
    data = load_two_parameter_sweep(f"autocat_abcd/barrido_k{i}_k{i+2}.pkl")
    
    r2_sweep = generate_r2_grid(data, k_change)
    with open(f"autocat_abcd/r2_sweep_k{i}_k{i+2}.pkl", "ab") as f:

        try:
            pickle.dump(r2_sweep, f)
        except Exception as e:
            print(f"An error has occured: {e}")

        finally: 
            print("Simulation run completed")

    data = None

        
    


