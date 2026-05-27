# Catalytic World

This repository contains the computational framework developed to simulate the chemical evolution and physical growth of volume-variable protocells. 

Specifically, this project investigates the kinetic requirements for the emergence and growth of early self-replicating chemical systems, focusing on the dynamics of a **minimal autocatalytic cycle**. The framework compares explicit multi-step chemical networks involving catalytic intermediates against widely used simplified third-order approximations, demonstrating the critical boundaries of these approximations in prebiotic chemistry.

## Features

* **Deterministic modeling (ODEs):** Simulates the macroscopic evolution of species concentrations and protocell volume using explicit multi-step reactions and lumped third-order kinetics.
* **Stochastic modeling (Gillespie algorithm):** Captures the inherent chemical noise of low-copy-number prebiotic networks to accurately model stochastic volume growth and division.

## Repository Structure


```text
catalytic_world/
│
├── src/                         # Core computational and visualization modules
│   ├── fun_gilles.py            # Stochastic simulation algorithms (Gillespie)
│   └── plotting.py              # Visualization
│
├── simulations/                 # Main execution scripts for experiments
│   ├── kinetic_models_comparison.py # Compares deterministic 3rd-order vs intermediate models
│   ├── ab0_sweep.py             # Catalyst initial concentration sweeps
│   ├── k_sweep.py               # 1D kinetic parameter sweeps
│   ├── k_sweep_2d.py            # 2D kinetic parameter sweeps (Alpha / R^2 heatmaps)
│   └── xyc_sweeps.py            # Sweeps for the XYC intermediate catalytic association
│
├── examples/                    # Reaction network definition files (.txt)
│
├── figures_TFM/                 # Generated output directory for thesis figures
├── figures_internship/          # Generated output directory for preliminary report figures
├── data_archive/                # Output directory for saved .pkl simulation results (Git ignored)
│
├── requirements.txt             # Python dependencies
└── README.md                    # Project documentation
