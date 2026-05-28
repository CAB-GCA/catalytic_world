import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import pandas as pd
import matplotlib.ticker as ticker
import os
import matplotlib.lines as mlines


def set_publication_style():
    """
    Sets global matplotlib parameters for thesis-quality figures.
    Call this once at the beginning of your main execution scripts.
    """
    plt.rcParams.update({
        'font.size': 14,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14,
        'lines.linewidth': 2.0,
        'figure.figsize': (8, 6),
        'figure.dpi': 300,
        'font.family': 'calibri',
        'axes.grid': False,
        'grid.alpha': 0.3,
        'image.cmap': 'viridis', 
        'mathtext.default':'regular'
    })


def plot_internship_abundances(times, abundances, species, fig_name=None):
    """
    Plots pure abundances (or concentrations if V=1) for the internship report.
    Automatically saves to 'figures_internship/'.
    """
    colors = plt.get_cmap('Set2').colors
    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    
    for i in range(len(species)):
        ax.plot(times, abundances[:, i], label=species[i], color=colors[i], linewidth=2, alpha=0.9)
        
    ax.set_ylabel("Concentration / Abundance")
    ax.set_xlabel("Time")
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    
    if fig_name is not None:
        import os
        os.makedirs("../figures_internship", exist_ok=True)
        plt.savefig(f"../figures_internship/{fig_name}.png", dpi=300)
    plt.show()


def plot_internship_conc_and_vol(times, abundances, volumes, species, fig_name=None):
    """
    Plots Concentration (Top, 1.5x height) and Volume (Bottom, 1x height).
    Automatically saves to 'figures_internship/'.
    """
    colors = plt.get_cmap('Set2').colors
    
    # 3:2 height ratio matches your [1.5, 1] request
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 5), 
                             gridspec_kw={'height_ratios': [1.5, 1]}, 
                             constrained_layout=True)
    
    ax1 = axes[0]
    for i in range(len(species)):
        ax1.plot(times, abundances[:, i] / volumes, label=species[i], color=colors[i], linewidth=1.6)
        
    ax1.set_ylabel("Concentration")
    ax1.set_xlabel("Time")
    ax1.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    
    ax2 = axes[1]
    ax2.plot(times, volumes, color='gray')
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Volume")
    
    if fig_name is not None:
        import os
        os.makedirs("../figures_internship", exist_ok=True)
        plt.savefig(f"../figures_internship/{fig_name}.png", dpi=300)
    plt.show()

def plot_time_evolution(times, abundances, volumes, species_indices, custom_labels, title="Time Evolution", color_map=None):
    """
    Plots the concentration of specific species over time.
    """
    fig, ax = plt.subplots()
    
    if color_map is None:
        colors = plt.cm.tab10(np.linspace(0, 1, len(species_indices)))
    else:
        colors = color_map(np.linspace(0, 1, len(species_indices)))
    
    for idx, (species_col, label) in enumerate(zip(species_indices, custom_labels)):
        concentration = abundances[:, species_col] / volumes
        ax.plot(times, concentration, label=label, color=colors[idx])
        
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Concentration")
    ax.set_title(title)
    
    # Legend is placed outside if there are many species
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
    fig.tight_layout()
    
    return fig, ax

def plot_sweep_with_replicates(param_values, raw_replicate_data, param_name, custom_label, log_x=True):
    """
    Plots the steady-state volume or concentration across a parameter sweep.
    """
    fig, ax = plt.subplots()
    
    # Calculate mean and standard deviation across replicates (axis=1)
    means = np.mean(raw_replicate_data, axis=1)
    stdevs = np.std(raw_replicate_data, axis=1)
    
    # Plotting with error bars
    ax.errorbar(param_values, means, yerr=stdevs, marker='o', linestyle='-', 
                color='teal', ecolor='gray', capsize=4, label=custom_label)
    
    if log_x:
        ax.set_xscale('log')
        ax.xaxis.set_minor_locator(ticker.NullLocator()) # REMOVE MINOR TICKS
        
    ax.set_xlabel(param_name)
    ax.set_ylabel("Steady State Value")
    ax.legend(loc='best')
    
    fig.tight_layout()
    return fig, ax

def plot_2d_sweep_heatmap(param_x, param_y, z_matrix, xlabel, ylabel, zlabel="Steady State Volume"):
    """
    Creates a heatmap for 2D parameter sweeps (e.g., k1 vs k3).
    """
    fig, ax = plt.subplots()
    
    # Create the heatmap
    cax = ax.pcolormesh(param_x, param_y, z_matrix, shading='auto')
    
    # Add a colorbar with its specific label
    cbar = fig.colorbar(cax, ax=ax)
    cbar.set_label(zlabel)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_minor_locator(ticker.NullLocator()) # REMOVE MINOR TICKS
    ax.yaxis.set_minor_locator(ticker.NullLocator()) # REMOVE MINOR TICKS

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    fig.tight_layout()
    return fig, ax

def save_figure(fig, filename):
    """
    Standardized saving function to ensure no figures are cut off.
    """
    fig.savefig(f"../figures_TFM/{filename}", bbox_inches='tight', dpi=300)
    plt.close(fig)
    
    
def _get_alpha_single_run(cond:float, result:tuple, values_used=None):
    """Helper function to calculate Alpha and R^2 for a single run."""
    _, times, volumes = result
    times = np.array(times)
    volumes = np.array(volumes)
    
    if values_used is None:
        start_idx = round(len(volumes) / 5)
        log_v = np.log(volumes[start_idx:])
        t_subset = times[start_idx:]
    else:
        values_used_idx = round(len(volumes) * values_used)
        log_v = np.log(volumes[:values_used_idx])
        t_subset = times[:values_used_idx]
    
    coefficients = np.polyfit(t_subset, log_v, 1)
    m_fit = coefficients[0]
    log_A_fit = coefficients[1]
    
    log_v_predicted = log_A_fit + m_fit * t_subset
    SS_res = np.sum((log_v - log_v_predicted)**2)
    SS_tot = np.sum((log_v - np.mean(log_v))**2)
    R_squared = 1 - (SS_res / SS_tot) if SS_tot != 0 else 0
    
    return {'Alpha': m_fit, 'R^2': R_squared}


def plot_combined_AB(results_list, xlim=None, ylim=None, lim=None, rep=0, fig_name=None):
    """
    results_list: List of 3 result dictionaries.
    Plots a 2x3 grid: Top row is Alpha/R^2 Sweep, Bottom row is Volume Evolution.
    """
    title_fs = 16
    label_fs = 15
    tick_fs = 14
    legend_fs = 14
    fig, axes = plt.subplots(2, 3, figsize=(12, 6), constrained_layout=True)
    
    titles_sweep = [r"(A)" +"\t"+r"$k_2 = 10^{-4}$", r"(B)"+"\t"+r"$k_2 = 1$", r"(C)"+"\t"+r"$k_2 = 10^4$"]
    titles_vol = [r"(D)"+"\t"+r"$k_2 = 10^{-4}$", r"(E)"+"\t"+r"$k_2 = 1$", r"(F)"+"\t"+r"$k_2 = 10^4$"]
    
    for i, results in enumerate(results_list):
        ax_sweep = axes[0, i]
        ax_vol = axes[1, i]

        # ==========================================
        # TOP ROW: BARRIDO / SWEEP (Alpha & R^2)
        # ==========================================
        res_copy = results.copy()
        if 0.0 in res_copy: res_copy.pop(0.0)
        
        stats_results = []
        for condition, replicates in res_copy.items():
            alphas, r2_values = [], []
            reps_list = replicates if isinstance(replicates, (list, tuple)) else [replicates]
            
            for r in reps_list:
                data_tuple = list(r.values())[0] if isinstance(r, dict) else r
                res = _get_alpha_single_run(condition, data_tuple, values_used=1) 
                alphas.append(res['Alpha'])
                r2_values.append(res['R^2'])
            
            stats_results.append({
                'Condition': condition,
                'Alpha_mean': np.mean(alphas),
                'Alpha_std': np.std(alphas),
                'R2_mean': np.mean(r2_values)
            })

        df = pd.DataFrame(stats_results).sort_values('Condition')

        # Primary Axis: Alpha
        ax_sweep.set_xlabel(r'$[\text{AB}]_0$', fontsize=label_fs) 
        ax_sweep.set_ylabel(r'$\alpha$', color="teal", fontsize=label_fs)
        ax_sweep.errorbar(df['Condition'], df['Alpha_mean'], yerr=df['Alpha_std'], 
                          color="teal", capsize=3, label=r'Mean $\alpha$')
        ax_sweep.tick_params(axis='y', labelcolor="teal", labelsize=tick_fs)
        ax_sweep.set_title(titles_sweep[i], fontsize=title_fs, loc='left')
        ax_sweep.axhline(y=0, color="teal", alpha=0.5, linestyle=":")
        
        
        ax_sweep.tick_params(axis='x', labelsize=tick_fs)
        # ax_sweep.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5)) # MaxNLocator can conflict with log scales occasionally
        ax_sweep.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))

        # Secondary Axis: R^2
        ax2_sweep = ax_sweep.twinx() 
        ax2_sweep.set_ylabel(r'$R^2$', color="orangered", fontsize=label_fs)
        ax2_sweep.scatter(df['Condition'], df['R2_mean'], rotation=270,labelpad=20,color="orangered", marker='s', alpha=0.6)
        ax2_sweep.tick_params(axis='y', labelcolor="orangered", labelsize=tick_fs)
        ax2_sweep.set_ylim((0, 1.05))
        ax2_sweep.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))

        # ==========================================
        # BOTTOM ROW: VOLUMES EVOLUTION
        # ==========================================
        if lim is not None:
            filtered_conditions = [c for c in results.keys() if lim[0] < float(c) < lim[1]]
        else:
            filtered_conditions = list(results.keys())
            
        filtered_conditions = sorted(filtered_conditions, key=float)
        
        if len(filtered_conditions) > 1:
            filtered_conditions.pop(0)
            filtered_conditions.pop(0)

        if filtered_conditions:
            color_gradient = plt.cm.rainbow(np.linspace(0, 1, len(filtered_conditions)))

            for j, condition in enumerate(filtered_conditions):
                if j % 2 != 0:  
                    continue
                replicates = results[condition]
                
                if not isinstance(replicates, tuple):
                    if isinstance(replicates[rep], dict):
                        data_tuple = list(replicates[rep].values())[0]
                    else:
                        data_tuple = replicates[rep]
                    a, t, v = data_tuple
                else:
                    a, t, v = replicates
                    
                ax_vol.plot(t, v, color=color_gradient[j], label=condition, alpha=0.9, linewidth=1.5)

        ax_vol.set_title(titles_vol[i], fontsize=title_fs, loc='left')
        ax_vol.set_yscale('log')
        ax_vol.yaxis.set_minor_locator(ticker.NullLocator()) # REMOVE MINOR TICKS
        
        ax_vol.set_xlabel("Time", fontsize=label_fs)
        ax_vol.set_ylabel("Volume", fontsize=label_fs)
        ax_vol.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
        ax_vol.tick_params(axis='both', labelsize=tick_fs)
        
        if i == 2:
            ax_vol.legend(title=rf"$[\text{{AB}}]_0$", fontsize=legend_fs-1, title_fontsize=title_fs-1,
                          loc='upper left', bbox_to_anchor=(1.05, 1.15))
            
        if xlim: ax_vol.set_xlim(0, xlim[i])
        if ylim: ax_vol.set_ylim(ylim)
            
    if fig_name:
        os.makedirs("../figures_TFM", exist_ok=True)
        plt.savefig(f"../figures_TFM/{fig_name}.png", dpi=300)
        
    plt.show()
    
def load_streamed_pickle(filename):
    data = {}
    try:
        with open(filename, "rb") as file:
            while True:
                try:
                    chunk = pickle.load(file)
                    data.update(chunk)
                except EOFError:
                    break
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    
    print(f"Successfully loaded {len(data)} simulation results from {filename}.")
    return data

def get_alpha_single_run(cond:float, result:tuple, values_used=None):
    _, times, volumes = result
    times = np.array(times)
    volumes = np.array(volumes)
    
    if values_used is None:
        start_idx = round(len(volumes) / 5)
        log_v = np.log(volumes[start_idx:])
        t_subset = times[start_idx:]
    else:
        values_used_idx = round(len(volumes) * values_used)
        log_v = np.log(volumes[:values_used_idx])
        t_subset = times[:values_used_idx]
    
    # Linear regression
    coefficients = np.polyfit(t_subset, log_v, 1)
    m_fit = coefficients[0]
    log_A_fit = coefficients[1]
    
    # R-squared calculation
    log_v_predicted = log_A_fit + m_fit * t_subset
    SS_res = np.sum((log_v - log_v_predicted)**2)
    SS_tot = np.sum((log_v - np.mean(log_v))**2)
    R_squared = 1 - (SS_res / SS_tot) if SS_tot != 0 else 0
    
    return {'Alpha': m_fit, 'R^2': R_squared}

def extract_stats_from_replicates(results_dict, values_used=1):
    res_copy = results_dict.copy()
    if 0.0 in res_copy: 
        res_copy.pop(0.0)
        
    stats_results = []
    
    for condition, replicates in res_copy.items():
        alphas, r2_values = [], []
        reps_list = replicates if isinstance(replicates, (list, tuple)) else [replicates]
        
        for r in reps_list:
            data_tuple = list(r.values())[0] if isinstance(r, dict) else r
            res = get_alpha_single_run(condition, data_tuple, values_used=values_used) 
            alphas.append(res['Alpha'])
            r2_values.append(res['R^2'])
            
        stats_results.append({
            'Condition': condition,
            'Alpha_mean': np.mean(alphas),
            'Alpha_std': np.std(alphas),
            'R2_mean': np.mean(r2_values)
        })
        
    return pd.DataFrame(stats_results).sort_values('Condition')



def plot_concentration_evolution(times, abundances, volumes, species, fig_name=None, xlim=None):
    """Plots the concentration evolution over time (Stochastic)."""
    colors = plt.get_cmap('Set2').colors
    fig, ax = plt.subplots(figsize=(5, 3), constrained_layout=True)
    
    for i in range(len(species)):
        ax.plot(times, abundances[:, i] / volumes, label=species[i].upper(), color=colors[i], alpha=0.9)
        
    ax.set_xlabel("Time")
    ax.set_ylabel("Concentration")
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), title="Species")
    
    if xlim is not None:
        ax.set_xlim((0, xlim))
        
    if fig_name is not None:
        os.makedirs("../figures_TFM", exist_ok=True)
        plt.savefig(f"../figures_TFM/{fig_name}.png", dpi=300)
    plt.show()

def plot_comparison_concentration(deterministic, stochastic, species, fig_name=None, xlim=None):
    """Plots deterministic vs stochastic concentration evolution."""
    colors = plt.get_cmap('Set2').colors
    con_det, t_det = deterministic
    a_sto, t_sto, v_sto = stochastic
    
    fig, ax = plt.subplots(figsize=(5, 2.5), constrained_layout=True)
    
    for i in range(len(species)):
        ax.plot(t_det, con_det[i], label=f"{species[i]} (Deterministic)", color=colors[i], alpha=0.9, linestyle='--')
        ax.plot(t_sto, a_sto[:, i] / v_sto, label=f"{species[i]} (Stochastic)", color=colors[i], alpha=0.5)
        
    ax.set_xlabel("Time")
    ax.set_ylabel("Concentration")

    # Custom Legend
    line_det = Line2D([0], [0], color='grey', linestyle='--', label='Deterministic')
    line_stoch = Line2D([0], [0], color='grey', linestyle='-', label='Stochastic')
    species_handles = [Patch(color=c, label=name) for c, name in zip(colors, species)]
    handles = species_handles + [Patch(color='none', label='')] + [line_det, line_stoch]

    ax.legend(handles=handles, loc="center right")
    
    if xlim is not None:
        ax.set_xlim((0, xlim))
    if fig_name is not None:
        os.makedirs("../figures_TFM", exist_ok=True)
        plt.savefig(f"../figures_TFM/{fig_name}.png", dpi=300)
    plt.show()

def plot_comparison_conc_and_vol(deterministic, stochastic, species, fig_name=None, xlim=None):
    """Plots deterministic vs stochastic for both Concentration (A) and Volume (B)."""
    colors = plt.get_cmap('Set2').colors
    con_det, t_det, v_det = deterministic
    a_sto, t_sto, v_sto = stochastic
    fig, axes = plt.subplots(figsize=(5, 5), nrows=2, ncols=1, constrained_layout=True)
    
    # Subplot A: Concentration
    ax1 = axes[0]
    for i in range(len(species)):
        ax1.plot(t_det, con_det[i], color=colors[i], alpha=0.9, linestyle='--')
        ax1.plot(t_sto, a_sto[:, i] / v_sto, color=colors[i], alpha=0.5)
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Concentration")
    ax1.set_title("(A)", loc="left", fontsize=16)

    # Subplot B: Volume
    ax2 = axes[1]
    ax2.plot(t_sto, v_sto, color='black', alpha=0.5, label='Volume (Stochastic)')
    ax2.plot(t_det, v_det, color='black', alpha=0.9, linestyle='--', label='Volume (Deterministic)')
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Volume")
    
    # REMOVE MINOR TICKS
    ax2.set_yscale('log')
    ax2.yaxis.set_minor_locator(ticker.NullLocator()) 

    ax2.set_title("(B)", loc="left", fontsize=16)

    # Custom Legend for Ax1
    line_det = Line2D([0], [0], color='grey', linestyle='--', label='Deterministic')
    line_stoch = Line2D([0], [0], color='grey', linestyle='-', label='Stochastic')
    species_handles = [Patch(color=c, label=name) for c, name in zip(colors, species)]
    handles = species_handles + [Patch(color='none', label='')] + [line_det, line_stoch]
    ax1.legend(handles=handles, loc='center right')
    
    if xlim is not None:
        ax1.set_xlim((0, xlim))
        ax2.set_xlim((0, xlim))
        
    if fig_name is not None:
        os.makedirs("../figures_TFM", exist_ok=True)
        plt.savefig(f"../figures_TFM/{fig_name}.png", dpi=300)
    plt.show()


def plot_precalculated_alpha(alpha_grid, k_values, k_in, k_out, min_log=-5, max_log=7, figname=None):
    """Plots the Alpha (Growth Rate) 2D Grid."""
    tick_fs = 14
    log_k = np.log10(k_values)
    mask = (log_k >= min_log) & (log_k <= max_log)
    
    floor_val = -6.0
    log_alpha_grid = np.full_like(alpha_grid, floor_val)
    np.log10(alpha_grid, where=(alpha_grid > 0), out=log_alpha_grid)
    
    plt.figure(figsize=(4, 3))
    filtered_k = log_k[mask]
    extent = [filtered_k[0], filtered_k[-1], filtered_k[0], filtered_k[-1]]

    current_cmap = plt.colormaps.get_cmap('plasma').copy()
    current_cmap.set_under('black') 
    current_cmap.set_bad('black') 
    
    vmin = floor_val + 0.1 
    vmax = np.max(log_alpha_grid) if np.max(log_alpha_grid) > vmin else vmin + 1
    
    im = plt.imshow(log_alpha_grid, origin='lower', extent=extent,
                    aspect='auto', cmap=current_cmap, vmin=vmin, vmax=vmax)

    fmt = ticker.FuncFormatter(lambda x, pos: f"$10^{{{int(x)}}}$")
    cbar = plt.colorbar(im, extend='min', format=fmt)
    
    desired_ticks = np.arange(np.ceil(vmin), np.floor(vmax) + 1, 2)
    all_ticks = np.sort(np.append(desired_ticks, floor_val))
    cbar.set_ticks(all_ticks)
    
    labels = ["$0$" if t == floor_val else f"$10^{{{int(t)}}}$" for t in all_ticks]
    cbar.set_ticklabels(labels)
    cbar.ax.tick_params(labelsize=tick_fs)
    cbar.set_label(r'Growth Rate ($\alpha$)', rotation=270, labelpad=20)
    
    ax = plt.gca()
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)
    
    plt.axvline(x=0, color='lightblue', linestyle='--', alpha=0.7)
    plt.xlabel(fr'$k_{{{k_in}}}$')
    plt.ylabel(fr'$k_{{{k_out}}}$')
    plt.tick_params(axis='both', labelsize=tick_fs)
    
    if figname:
        os.makedirs("../figures_TFM", exist_ok=True)
        plt.savefig(f"../figures_TFM/{figname}.png", bbox_inches='tight', dpi=300)
    plt.tight_layout()
    plt.show()


def plot_r2_map(r2_grid, k_values, k_in, k_out, min_log=-5, max_log=7, figname=None):
    """Plots the pre-calculated R^2 grid as a heatmap."""
    tick_fs = 14
    log_k = np.log10(k_values)
    mask = (log_k >= min_log) & (log_k <= max_log)
    
    filtered_grid = r2_grid[mask][:, mask]
    filtered_k = log_k[mask]
    
    # --- THE FIX ---
    # We only check if the R2 grid is empty here!
    if filtered_grid.size == 0:
        print("No data in range.")
        return
        
    plt.figure(figsize=(4, 3))
    extent = [filtered_k[0], filtered_k[-1], filtered_k[0], filtered_k[-1]]
    
    im = plt.imshow(filtered_grid, origin='lower', extent=extent, 
                    aspect='auto', cmap='viridis', vmin=0, vmax=1)
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize=tick_fs)
    cbar.set_label(r'$R^2$ of the fit', rotation=270, labelpad=20)
    
    fmt = ticker.FuncFormatter(lambda x, pos: f"$10^{{{int(x)}}}$")
    ax = plt.gca()
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)
    
    plt.xlabel(fr'$k_{{{k_in}}}$')
    plt.ylabel(fr'$k_{{{k_out}}}$')
    plt.tick_params(axis='both', labelsize=tick_fs)
    
    if figname:
        os.makedirs("../figures_TFM", exist_ok=True)
        plt.savefig(f"../figures_TFM/{figname}.png", bbox_inches='tight', dpi=300)
    plt.tight_layout()
    plt.show()


    
def plot_both_2d_sweeps(alpha_grid, r2_grid, k_values, k_in, k_out, min_log=-np.inf, max_log=np.inf, figname=None):
    """Plots both the Alpha (A) and R^2 (B) heatmaps side-by-side."""
    import matplotlib.ticker as ticker  # <--- THE FIX: Explicitly load ticker here

    tick_fs = 14
    title_fs = 16
    if r2_grid.size == 0 or alpha_grid.size == 0:
        print("No data to plot.")
        return

    log_k = np.log10(k_values)
    mask = (log_k >= min_log) & (log_k <= max_log)
    filtered_k = log_k[mask]
    
    if len(filtered_k) == 0:
        filtered_k = log_k 
        
    extent = [filtered_k[0], filtered_k[-1], filtered_k[0], filtered_k[-1]]
    fmt = ticker.FuncFormatter(lambda x, pos: f"$10^{{{int(x)}}}$")

    # --- 3-Column Layout ---
    fig, axes = plt.subplots(1, 3, figsize=(8,2.914), 
                             gridspec_kw={'width_ratios': [1, 1, 0.35]}, 
                             constrained_layout=True)
    ax1, ax2, ax3 = axes
    ax3.axis('off') # Hide the empty column

    # --- LEFT COLUMN: Alpha Grid ---
    floor_val = -6.0
    log_alpha_grid = np.full_like(alpha_grid, floor_val)
    np.log10(alpha_grid, where=(alpha_grid > 0), out=log_alpha_grid)
    
    current_cmap = plt.colormaps.get_cmap('plasma').copy()
    current_cmap.set_under('black') 
    current_cmap.set_bad('black') 
    
    vmin = floor_val + 0.1 
    vmax = np.max(log_alpha_grid) if np.max(log_alpha_grid) > vmin else vmin + 1
    
    im1 = ax1.imshow(log_alpha_grid, origin='lower', extent=extent,
                     aspect='auto', cmap=current_cmap, vmin=vmin, vmax=vmax)

    cbar1 = fig.colorbar(im1, ax=ax1, extend='min', format=fmt)
    desired_ticks = np.arange(np.ceil(vmin), np.floor(vmax) + 1, 2)
    all_ticks = np.sort(np.append(desired_ticks, floor_val))
    cbar1.set_ticks(all_ticks)
    labels = ["$0$" if t == floor_val else f"$10^{{{int(t)}}}$" for t in all_ticks]
    cbar1.set_ticklabels(labels)
    cbar1.ax.tick_params(labelsize=tick_fs)
    cbar1.set_label(r'Growth rate $\alpha$', rotation=270, labelpad=20)
    
    ax1.xaxis.set_major_formatter(fmt)
    ax1.yaxis.set_major_formatter(fmt)
    ax1.set_xlabel(fr'$k_{{{k_in}}}$')
    ax1.set_ylabel(fr'$k_{{{k_out}}}$')
    ax1.tick_params(axis='both', labelsize=tick_fs)
    ax1.set_title("(A)", loc="left", x=-0.05, y=1, pad=10, fontsize=title_fs)
    ax1.grid(False)

    # --- RIGHT COLUMN: R^2 Grid ---
    im2 = ax2.imshow(r2_grid, origin='lower', extent=extent, 
                     aspect='auto', cmap='viridis', vmin=0, vmax=1)
    
    cbar2 = fig.colorbar(im2, ax=ax2)
    cbar2.ax.tick_params(labelsize=tick_fs)
    cbar2.set_label(r'$R^2$ of the fit', rotation=270, labelpad=20)
    
    ax2.xaxis.set_major_formatter(fmt)
    ax2.yaxis.set_major_formatter(fmt)
    ax2.set_xlabel(fr'$k_{{{k_in}}}$')
    ax2.set_ylabel(fr'$k_{{{k_out}}}$')
    ax2.tick_params(axis='both', labelsize=tick_fs)
    ax2.set_title("(B)", loc="left", x=-0.05, y=1, pad=10, fontsize=title_fs)

    if figname:
        import os
        os.makedirs("../figures_TFM", exist_ok=True)
        plt.savefig(f"../figures_TFM/{figname}.png", bbox_inches='tight', dpi=300)
    plt.show()
    
def plot_combined_k_sweep_and_volumes(results: dict, k_i: int, xlim=None, ylim=None, lim=None, rep=0, figname=None):
    """
    Plots a 1x2 grid for a 1D kinetic parameter sweep.
    Left: Alpha and R^2 vs k_i. Right: Volume evolution over time.
    """
    import pandas as pd
    import matplotlib.ticker as ticker
    
    title_fs = 16
    label_fs = 15
    tick_fs = 14
    legend_fs = 14

    # --- THE FIX: 3-Column Layout ---
    # Width ratios of 1:1 force the first two plots to be perfectly equal.
    fig, axes = plt.subplots(1, 3, figsize=(8,2.914), 
                             gridspec_kw={'width_ratios': [1, 1, 0.35]}, 
                             constrained_layout=True)    
    ax_sweep = axes[0]
    ax_vol = axes[1]
    ax_leg = axes[2]
    ax_leg.axis('off') # Hide the third axis completely
    
    # ==========================================
    # LEFT COLUMN: Alpha & R^2 Sweep
    # ==========================================
    stats_results = []
    
    for condition, replicates in results.items():
        alphas = []
        r2_values = []
        
        reps_list = replicates if isinstance(replicates, (list, tuple)) else [replicates]
        
        for r in reps_list:
            data_tuple = list(r.values())[0] if isinstance(r, dict) else r
            res = _get_alpha_single_run(condition, data_tuple) 
            alphas.append(res['Alpha'])
            r2_values.append(res['R^2'])
        
        stats_results.append({
            'Condition': float(condition),
            'Alpha_mean': np.mean(alphas),
            'Alpha_std': np.std(alphas),
            'R2_mean': np.mean(r2_values)
        })

    df = pd.DataFrame(stats_results).sort_values('Condition')

    # Primary Axis: Alpha
    ax_sweep.set_xlabel(fr'$k_{{{k_i+1}}}$', fontsize=label_fs)
    ax_sweep.set_ylabel(r'$\alpha$', color="teal", fontsize=label_fs)
    ax_sweep.errorbar(df['Condition'], df['Alpha_mean'], yerr=df['Alpha_std'], 
                      color="teal", capsize=4, label=r'Mean $\alpha$')
    
    ax_sweep.tick_params(axis='y', labelcolor="teal", labelsize=tick_fs)
    ax_sweep.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
    ax_sweep.axhline(y=0, color="teal", alpha=0.6, linestyle=":")
    
    ax_sweep.set_xscale('log')
    ax_sweep.xaxis.set_minor_locator(ticker.NullLocator())
    
    ax_sweep.tick_params(axis='x', labelsize=tick_fs)
    ax_sweep.set_title("(C)", loc="left", x=-0.05, y=1, pad=10, fontsize=title_fs)
    ax_sweep.grid(False)

    # Secondary Axis: R^2
    ax2_sweep = ax_sweep.twinx()  
    ax2_sweep.set_ylabel(r'$R^2$', color="orangered", fontsize=label_fs)  
    ax2_sweep.scatter(df['Condition'], df['R2_mean'], color="orangered", marker='s', alpha=0.6)
    ax2_sweep.tick_params(axis='y', labelcolor="orangered", labelsize=tick_fs)
    ax2_sweep.set_ylim((0, 1.05))
    ax2_sweep.grid(False)
    
    # ==========================================
    # RIGHT COLUMN: Volumes vs Time
    # ==========================================
    if lim is not None:
        filtered_conditions = [c for c in results.keys() if lim[0] < float(c) < lim[1]] 
    else:
        filtered_conditions = list(results.keys())
        
    filtered_conditions = sorted(filtered_conditions, key=float)
    color_gradient = plt.cm.rainbow(np.linspace(0, 1, len(filtered_conditions)))

    for j, condition in enumerate(filtered_conditions):
        replicates = results[condition]
        formatted_label = f"${float(condition):.1e}$"
        
        if not isinstance(replicates, tuple):
            if isinstance(replicates[rep], dict):
                data_tuple = list(replicates[rep].values())[0]
            else:
                data_tuple = replicates[rep]
            a, t, v = data_tuple
        else:
            a, t, v = replicates
            
        ax_vol.plot(t, v, color=color_gradient[j], label=formatted_label, alpha=0.9, linewidth=2)    

    ax_vol.set_yscale('log')
    ax_vol.yaxis.set_minor_locator(ticker.NullLocator())
    
    ax_vol.set_xlabel("Time", fontsize=label_fs)
    ax_vol.set_ylabel("Volume", fontsize=label_fs)
    ax_vol.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
    ax_vol.tick_params(axis='both', labelsize=tick_fs)
    ax_vol.set_title("(D)", loc="left", x=-0.05, y=1, pad=10, fontsize=title_fs)

    if xlim is not None:
        ax_vol.set_xlim((0, xlim))
    if ylim is not None:
        ax_vol.set_ylim((5e1, ylim))

    # ==========================================
    # LEGEND INJECTION
    # ==========================================
    handles, labels = ax_vol.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    
    sorted_labels = sorted(by_label.keys(), key=lambda x: float(x.replace('$', '').replace('{', '').replace('}', '')))
    sorted_handles = [by_label[k] for k in sorted_labels]
    
    # Attach the legend to the invisible 3rd axis
    ax_leg.legend(sorted_handles, sorted_labels, title=rf"$k_{{{k_i+1}}}$", fontsize=legend_fs-2,
                  loc='center left', title_fontsize=title_fs)

    if figname:
        import os
        os.makedirs("../figures_TFM", exist_ok=True)
        plt.savefig(f"../figures_TFM/{figname}.png", dpi=300)
        
    plt.show()
    
    
def plot_master_4panel_figure(alpha_grid, r2_grid, k_values_2d, k_in, k_out, 
                              results_1d, k_i, xlim_1d=None, ylim_1d=None, lim_1d=None, 
                              min_log=-np.inf, max_log=np.inf, figname=None):
    """
    Creates a 2x2 master figure. 
    Top row: 2D Heatmaps (Alpha and R^2). 
    Bottom row: 1D parameter sweeps (Alpha/R^2 and Volume Evolution).
    """
    import pandas as pd
    import matplotlib.ticker as ticker
    from matplotlib.lines import Line2D
    import matplotlib.pyplot as plt

    # Global Font Sizes
    title_fs = 16
    label_fs = 15
    tick_fs = 14
    legend_fs = 14

    # --- SETUP THE 2x3 GRID ---
    fig, axes = plt.subplots(2, 3, figsize=(10,6.25), 
                             gridspec_kw={'width_ratios': [1, 1, 0.35], 'height_ratios': [1, 1]}, 
                             constrained_layout=True)
    
    ax_alpha_2d, ax_r2_2d, ax_empty_top = axes[0]
    ax_alpha_1d, ax_vol_1d, ax_leg_bot = axes[1]
    
    ax_empty_top.axis('off') # Spacer for top row
    ax_leg_bot.axis('off')   # Holder for bottom row legend

    fmt = ticker.FuncFormatter(lambda x, pos: f"$10^{{{int(x)}}}$")

    # ==========================================
    # TOP ROW: 2D HEATMAPS
    # ==========================================
    if r2_grid.size > 0 and alpha_grid.size > 0:
        log_k = np.log10(k_values_2d)
        mask = (log_k >= min_log) & (log_k <= max_log)
        filtered_k = log_k[mask] if len(log_k[mask]) > 0 else log_k
        extent = [filtered_k[0], filtered_k[-1], filtered_k[0], filtered_k[-1]]

        # --- (A) 2D Alpha Grid ---
        floor_val = -6.0
        log_alpha_grid = np.full_like(alpha_grid, floor_val)
        np.log10(alpha_grid, where=(alpha_grid > 0), out=log_alpha_grid)
        
        cmap_plasma = plt.colormaps.get_cmap('plasma').copy()
        cmap_plasma.set_under('black') 
        cmap_plasma.set_bad('black') 
        
        vmin = floor_val + 0.1 
        vmax = np.max(log_alpha_grid) if np.max(log_alpha_grid) > vmin else vmin + 1
        
        im1 = ax_alpha_2d.imshow(log_alpha_grid, origin='lower', extent=extent,
                                 aspect='auto', cmap=cmap_plasma, vmin=vmin, vmax=vmax)

        # THE FIX: Create a tightly bound box for the colorbar (x0, y0, width, height)
        cax1 = ax_alpha_2d.inset_axes([1.05, 0, 0.05, 1])
        cbar1 = fig.colorbar(im1, cax=cax1, extend='min', format=fmt)
        
        desired_ticks = np.arange(np.ceil(vmin), np.floor(vmax) + 1, 2)
        all_ticks = np.sort(np.append(desired_ticks, floor_val))
        cbar1.set_ticks(all_ticks)
        cbar1.set_ticklabels(["$0$" if t == floor_val else f"$10^{{{int(t)}}}$" for t in all_ticks])
        cbar1.ax.tick_params(labelsize=tick_fs)
        cbar1.set_label(r'Growth rate $\alpha$', rotation=270, labelpad=20)
        
        ax_alpha_2d.xaxis.set_major_formatter(fmt)
        ax_alpha_2d.yaxis.set_major_formatter(fmt)
        ax_alpha_2d.set_xlabel(fr'$k_{{{k_in}}}$', fontsize=label_fs)
        ax_alpha_2d.set_ylabel(fr'$k_{{{k_out}}}$', fontsize=label_fs)
        ax_alpha_2d.tick_params(axis='both', labelsize=tick_fs)
        ax_alpha_2d.set_title("(A)", loc="left", x=-0.05, y=1, pad=10, fontsize=title_fs)

        # --- (B) 2D R^2 Grid ---
        im2 = ax_r2_2d.imshow(r2_grid, origin='lower', extent=extent, 
                              aspect='auto', cmap='viridis', vmin=0, vmax=1)
        
        # THE FIX: Create a tightly bound box for the colorbar
        cax2 = ax_r2_2d.inset_axes([1.05, 0, 0.05, 1])
        cbar2 = fig.colorbar(im2, cax=cax2)
        
        cbar2.ax.tick_params(labelsize=tick_fs)
        cbar2.set_label(r'$R^2$', rotation=270, labelpad=20)
        
        ax_r2_2d.xaxis.set_major_formatter(fmt)
        ax_r2_2d.yaxis.set_major_formatter(fmt)
        ax_r2_2d.set_xlabel(fr'$k_{{{k_in}}}$', fontsize=label_fs)
        ax_r2_2d.set_ylabel(fr'$k_{{{k_out}}}$', fontsize=label_fs)
        ax_r2_2d.tick_params(axis='both', labelsize=tick_fs)
        ax_r2_2d.set_title("(B)", loc="left", x=-0.05, y=1, pad=10, fontsize=title_fs)

    # ==========================================
    # BOTTOM ROW: 1D SWEEPS
    # ==========================================
    if results_1d:
        stats_results = []
        for condition, replicates in results_1d.items():
            alphas, r2_values = [], []
            reps_list = replicates if isinstance(replicates, (list, tuple)) else [replicates]
            
            for r in reps_list:
                data_tuple = list(r.values())[0] if isinstance(r, dict) else r
                res = _get_alpha_single_run(condition, data_tuple) 
                alphas.append(res['Alpha'])
                r2_values.append(res['R^2'])
            
            stats_results.append({
                'Condition': float(condition),
                'Alpha_mean': np.mean(alphas),
                'Alpha_std': np.std(alphas),
                'R2_mean': np.mean(r2_values)
            })

        df = pd.DataFrame(stats_results).sort_values('Condition')

        # --- (C) 1D Alpha & R^2 Sweep ---
        ax_alpha_1d.set_xlabel(fr'$k_{{{k_i+1}}}$', fontsize=label_fs)
        ax_alpha_1d.set_ylabel(r'Growth rate $\alpha$', color="teal", fontsize=label_fs)
        ax_alpha_1d.errorbar(df['Condition'], df['Alpha_mean'], yerr=df['Alpha_std'], 
                             color="teal", capsize=4, label=r'Mean $\alpha$')
        
        ax_alpha_1d.tick_params(axis='y', labelcolor="teal", labelsize=tick_fs)
        ax_alpha_1d.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
        ax_alpha_1d.axhline(y=0, color="teal", alpha=0.6, linestyle=":")
        
        ax_alpha_1d.set_xscale('log')
        ax_alpha_1d.xaxis.set_minor_locator(ticker.NullLocator())
        ax_alpha_1d.tick_params(axis='x', labelsize=tick_fs)
        ax_alpha_1d.set_title("(C)", loc="left", x=-0.05, y=1, pad=10, fontsize=title_fs)

        ax2_sweep = ax_alpha_1d.twinx()  
        ax2_sweep.set_ylabel(r'$R^2$', color="orangered",rotation=270,labelpad=20, fontsize=label_fs)  
        ax2_sweep.scatter(df['Condition'], df['R2_mean'], color="orangered", marker='s', alpha=0.6)
        ax2_sweep.tick_params(axis='y', labelcolor="orangered", labelsize=tick_fs)
        ax2_sweep.set_ylim((0, 1.05))
        
        # --- (D) 1D Volume Evolution ---
        if lim_1d is not None:
            filtered_conds = [c for c in results_1d.keys() if lim_1d[0] < float(c) < lim_1d[1]] 
        else:
            filtered_conds = list(results_1d.keys())
            
        filtered_conds = sorted(filtered_conds, key=float)
        color_grad = plt.cm.rainbow(np.linspace(0, 1, len(filtered_conds)))

        for j, condition in enumerate(filtered_conds):
            replicates = results_1d[condition]
            formatted_label = f"${float(condition):.1e}$"
            data_tuple = replicates[0] if isinstance(replicates, (list, tuple)) else replicates
            if isinstance(data_tuple, dict): data_tuple = list(data_tuple.values())[0]
            _, t, v = data_tuple
                
            ax_vol_1d.plot(t, v, color=color_grad[j], label=formatted_label, alpha=0.9, linewidth=2)    

        ax_vol_1d.set_yscale('log')
        ax_vol_1d.yaxis.set_minor_locator(ticker.NullLocator())
        ax_vol_1d.set_xlabel("Time", fontsize=label_fs)
        ax_vol_1d.set_ylabel("Volume", fontsize=label_fs)
        ax_vol_1d.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
        ax_vol_1d.tick_params(axis='both', labelsize=tick_fs)
        ax_vol_1d.set_title("(D)", loc="left", x=-0.05, y=1, pad=10, fontsize=title_fs)

        if xlim_1d is not None: ax_vol_1d.set_xlim((0, xlim_1d))
        if ylim_1d is not None: ax_vol_1d.set_ylim((5e1, ylim_1d))

        # --- Inject Legend into ax_leg_bot ---
        handles, labels = ax_vol_1d.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        sorted_labels = sorted(by_label.keys(), key=lambda x: float(x.replace('$', '').replace('{', '').replace('}', '')))
        sorted_handles = [by_label[k] for k in sorted_labels]
        
        ax_leg_bot.legend(sorted_handles, sorted_labels, title=rf"$k_{{{k_i+1}}}$", 
                          fontsize=legend_fs-2, loc='center left', title_fontsize=title_fs)

    if figname:
        import os
        os.makedirs("../figures_TFM", exist_ok=True)
        plt.savefig(f"../figures_TFM/{figname}.png", bbox_inches='tight', dpi=300)
        
    plt.show()

def plot_deterministic_dual_comparison(t_to, conc_to, v_to, species_to, 
                                       t_int_eq, conc_int_eq, v_int_eq, 
                                       t_int_br, conc_int_br, v_int_br, species_int, figname=None):
    """Plots a 2x2 grid. Left: Equivalent parameters. Right: Broken assumptions."""
    
    cmap_colors = list(plt.get_cmap('Set2').colors)
    linestyles = ["-", "--"]

    fig, axes = plt.subplots(figsize=(8.5,4), nrows=2, ncols=2, height_ratios=[2,1], 
                             constrained_layout=True) 
    
    species_color_map = {}
    
    # helper to plot one column
    def plot_column(ax_conc, ax_vol, t_int, conc_int, v_int, col_titles):
        # Concentration
        for i in range(len(species_to)):
            color = cmap_colors[i % len(cmap_colors)]
            species_color_map[species_to[i]] = color
            ax_conc.plot(t_to, conc_to[i], color=color, alpha=0.5, linestyle=linestyles[0], linewidth=2)
                     
        for i in range(len(species_int)):
            color = cmap_colors[i % len(cmap_colors)]
            species_color_map[species_int[i]] = color
            ax_conc.plot(t_int, conc_int[i], color=color, alpha=0.7, linestyle=linestyles[1], linewidth=2.5)
                     
        ax_conc.set_ylabel("Concentration")
        ax_conc.set_xlabel("Time")
        ax_conc.set_title(col_titles[0], loc="left")

        # Volume
        ax_vol.plot(t_to, v_to, color='black', linestyle=linestyles[0], alpha=0.5, linewidth=2)
        ax_vol.plot(t_int, v_int, color='black', linestyle=linestyles[1], alpha=0.7, linewidth=2.5)
        ax_vol.set_xlabel("Time")
        ax_vol.set_ylabel("Volume")
        ax_vol.set_title(col_titles[1], loc="left")
        
        ax_vol.set_yscale('log')
        ax_vol.yaxis.set_minor_locator(ticker.NullLocator()) # REMOVE MINOR TICKS

    # Plot Left Column (Equivalent)
    plot_column(axes[0, 0], axes[1, 0], t_int_eq, conc_int_eq, v_int_eq, ["(A)","(C)"])
    # Plot Right Column (Broken)
    plot_column(axes[0, 1], axes[1, 1], t_int_br, conc_int_br, v_int_br, ["(B)","(D)"])
    
    # Custom Legend (Shared for entire figure)
    species_handles = [Patch(color=color, label=sp) for sp, color in species_color_map.items()]
    species_handles.append(Patch(color='black', label='Volume'))
    line_to = Line2D([0], [0], color='gray', linestyle=linestyles[0], alpha=0.5, linewidth=2, label='Third order')
    line_int = Line2D([0], [0], color='gray', linestyle=linestyles[1], alpha=0.7, linewidth=2.5, label='Intermediates')
    
    handles = species_handles + [Patch(color='none', label='')] + [line_to, line_int]
    axes[0, 1].legend(handles=handles, loc='upper left', bbox_to_anchor=(1.05, 1.05))
    
    if figname:
        os.makedirs("../figures_TFM", exist_ok=True)
        plt.savefig(f"../figures_TFM/{figname}.png", dpi=300)
    plt.show()


def plot_volume_sweep_dual_comparison(results_to, results_int_eq, results_int_br, k_values, figname=None):

    
    colors_sweep = plt.cm.rainbow(np.linspace(0, 1, len(k_values)))
    fig, axes = plt.subplots(figsize=(8.5,3.5), ncols=2, constrained_layout=True)

    def plot_sweep_axis(ax, res_int, col_title):
        for n, k_val in enumerate(k_values):
            t_to, v_to = results_to[n]
            t_int, v_int = res_int[n]
            ax.plot(t_to, v_to, label=fr"$k_1$ = {k_val:.1e} (third-order)", 
                    color=colors_sweep[n], alpha=0.4, linewidth=2)
            ax.plot(t_int, v_int, label=fr"$k_3$ = {k_val:.1e} (intermediate)", 
                    color=colors_sweep[n], alpha=0.7, linewidth=2.5, linestyle='--')
        ax.set_xlabel("Time")
        ax.set_ylabel("Volume")
        
        ax.set_yscale('log')
        ax.yaxis.set_minor_locator(ticker.NullLocator()) # REMOVE MINOR TICKS
        
        ax.set_title(col_title, loc="left")

    plot_sweep_axis(axes[0], results_int_eq, "(E)")
    plot_sweep_axis(axes[1], results_int_br, "(F)")

    # Custom Legend
    line_to = mlines.Line2D([], [], color='gray', linestyle='-', label='Third-order')
    line_int = Line2D([], [], color='gray', linestyle='--', label='Intermediates')
    species_handles = [Patch(color=c, label=f'${val:.1e}$') for c, val in zip(colors_sweep, k_values)]
    handles = species_handles + [Patch(color='none', label='')] + [line_to, line_int]

    axes[1].legend(handles=handles, loc='center left', bbox_to_anchor=(1.05, 0.5), 
                   title=r"$k_1^{T}$ and $k_{eff}^{Int}$")

    if figname:
        os.makedirs("../figures_TFM", exist_ok=True)
        plt.savefig(f"../figures_TFM/{figname}.png", dpi=300)
    plt.show()