# Catalytic World

This repository contains the computational framework developed to simulate the chemical evolution and physical growth of volume-variable protocells. 

Specifically, this project investigates the kinetic requirements for the emergence and stability of early self-replicating chemical systems, focusing on the dynamics of a **minimal autocatalytic cycle**. The framework compares explicit multi-step chemical networks involving catalytic intermediates against widely used simplified third-order approximations, demonstrating the critical boundaries of these approximations in prebiotic chemistry.

## Features

* **Deterministic modeling (ODEs):** Simulates the macroscopic evolution of species concentrations and protocell volume using explicit multi-step reactions and lumped third-order kinetics.
* **Stochastic modeling (Gillespie algorithm):** Captures the inherent chemical noise of low-copy-number prebiotic networks to accurately model stochastic volume growth and division.

## Repository Structure

```text
catalytic_world/
│
├── src/
│   ├── plotting.py              # Visualization suite for time evolution and parameter sweeps
│   └── fun_gilles.py            # Stochastic simulation algorithms (Gillespie) 
│
├── kinetic_equivalence.py       # Main script to compare deterministic third-order vs. intermediate models
│
├── figures_TFM/                 # Generated output directory for thesis figures
├── figures_internship/          # Generated output directory for preliminary report figures
│
├── requirements.txt             # Python dependencies
└── README.md                    # Project documentation
