import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')  # allow importing from parent directory
from fun_gilles import *
import pickle
import pandas as pd

colors = plt.get_cmap('Set2').colors

def load_streamed_pickle(filename):
    data = {}
    try:
        with open(filename, "rb") as file:
            while True:
                try:
                    # Load one object (which is a dict of {k_i: (abundances, times, volumes)})
                    chunk = pickle.load(file)
                    # Merge the loaded chunk into the main dictionary
                    data.update(chunk)
                except EOFError:
                    # Reached the end of the file
                    break
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    
    print(f"Successfully loaded {len(data)} simulation results.")
    return data

def get_alpha_single_run(cond:float, result:tuple):
    _, times, volumes = result
        
    times = np.array(times)
    volumes = np.array(volumes)
    
    # Use data points from index len/5 onwards
    start_idx = round(len(volumes) / 5)
    log_v = np.log(volumes[start_idx:])
    t_subset = times[start_idx:]
    
    # 2. Perform linear regression
    # polyfit(x, y, 1) returns [slope, intercept]
    coefficients = np.polyfit(t_subset, log_v, 1)
    m_fit = coefficients[0]      # This is alpha (slope)
    log_A_fit = coefficients[1]  # This is the intercept
    A_fit = np.exp(log_A_fit)
    
    # --- R-SQUARED CALCULATION ---
    log_v_predicted = log_A_fit + m_fit * t_subset
    log_v_mean = np.mean(log_v)
    SS_res = np.sum((log_v - log_v_predicted)**2)
    SS_tot = np.sum((log_v - log_v_mean)**2)
    R_squared = 1 - (SS_res / SS_tot) if SS_tot != 0 else 0
    
    # Return a dictionary directly
    return {
        'Condition': cond,
        'Alpha': m_fit,
        'Scaling Const. (A)': A_fit,
        'R^2': R_squared
    }

results = load_streamed_pickle("drago_results/barrido_ab0_k1_0.0001.pkl")


# 1. Processing the statistical data
# Assuming barrido_ab_k1_em4 is a dict: {concentration: [(ab, t, v), (ab, t, v), ...]}
stats_results = []

for condition, replicates in results.items():
    alphas = []
    r2_values = []
    if not isinstance(replicates, tuple):
        for rep in replicates:
            # If rep is a dictionary {cond: (ab, t, v)}, get the tuple:
            if isinstance(rep, dict):
                # This gets the first value in the dict, which is your (ab, t, v) tuple
                data_tuple = list(rep.values())[0]
            else:
                data_tuple = rep
                
            res = get_alpha_single_run(condition, data_tuple) 
            alphas.append(res['Alpha'])
            r2_values.append(res['R^2'])
    else:
        res = get_alpha_single_run(condition, replicates) 
        alphas.append(res['Alpha'])
        r2_values.append(res['R^2'])
    
    stats_results.append({
        'Condition': condition,
        'Alpha_mean': np.mean(alphas),
        'Alpha_std': np.std(alphas),
        'R2_mean': np.mean(r2_values)
    })

# Convert to DataFrame or Sortable List for plotting
stats_results = pd.DataFrame(stats_results).sort_values('Condition')
print(stats_results)

# 2. Plotting with Error Bars
fig, ax1 = plt.subplots(figsize=(8, 5))

# Primary Axis: Alpha with Error Bars
ax1.set_xlabel(r'$[\text{AB}]_0$') 
ax1.set_ylabel(r'$\alpha$', color=colors[1])

ax1.errorbar(stats_results['Condition'], stats_results['Alpha_mean'], 
             yerr=stats_results['Alpha_std'], 
             fmt='o-', color=colors[1], capsize=4, label=r'Mean $\alpha$')

ax1.tick_params(axis='y', labelcolor=colors[1])
ax1.grid(True, linestyle=':', alpha=0.5, which='both')
ax1.axhline(y=0, color=colors[1], alpha=0.6, linestyle=":")
ax1.set_xscale('log')

# Secondary Axis: R^2
ax2 = ax1.twinx()  
ax2.set_ylabel(r'$R^2$', color=colors[2])  
ax2.scatter(stats_results['Condition'], stats_results['R2_mean'], 
            color=colors[2], marker='s', alpha=0.6)

ax2.tick_params(axis='y', labelcolor=colors[2])
ax2.set_ylim((0, 1.05))
ax2.set_xscale('log')

plt.title(r'Growth rate $\alpha$ vs Initial Catalyst Concentration')
fig.tight_layout()
plt.show()