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

def plot_division_2d_sweep(raw_data, k_values, n_divisions=1):
    """
    results: The loaded pickle data (dictionary)
    k_values: The np.logspace array used for the sweep
    k1_idx, k2_idx: The indices of the constants being swept (e.g., 0 and 2)
    """
    n = len(k_values)
    # Initialize a grid to store the mean division times
    td_grid = np.full((n, n), np.nan)
      
    n = len(k_values)
    for entry in raw_data:
        for condition, replicates in entry.items():
            # Map flat 'condition' index back to row/col
            row = condition // n
            col = condition % n
            
            doubling_times = []
            for rep in replicates:
                data_tuple = list(rep.values())[0] if isinstance(rep, dict) else rep
                _, times, volumes = data_tuple
                
                v_initial = volumes[0]
                v_target = v_initial * (2 ** n_divisions)
                
                # Find the first time the volume crosses the threshold
                idx_doubled = np.where(volumes >= v_target)[0]
                
                if len(idx_doubled) > 0:
                    doubling_times.append(times[idx_doubled[0]])
                else:
                    doubling_times.append(np.nan)
            
            # Store the mean time in the grid
            if not np.all(np.isnan(doubling_times)):
                td_grid[row, col] = np.nanmean(doubling_times)
    
    return td_grid

k_change = np.logspace(-5,7,25)
n_div= 2

for i in [0, 1]:
    data = load_two_parameter_sweep(f"autocat_abcd/barrido_k{i}_k{i+2}.pkl")
    
    alpha_sweep = plot_division_2d_sweep(data, k_change, n_div)
    with open(f"autocat_abcd/division{n_div}_sweep_k{i}_k{i+2}.pkl", "ab") as f:

        try:
            pickle.dump(alpha_sweep, f)
        except Exception as e:
            print(f"An error has occured: {e}")

        finally: 
            print("Simulation run completed")
    
    data = None
