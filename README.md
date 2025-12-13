# Project-paper
Particle orbital simulation using numerical solvers and perturbation expressions

## Purpose
This project simulates motion of one particle under the influence og gravitational pull from the Sun,
pressure radiation force (radial component) and mass loss due to sputtering. Simulations are done
with Runge Kutta 4th order scheme by using scipy.integrate.solve_ivp and Leapfrog integrators as 
numerical solvers, and analytically by perturbed expressions for radial distance r, angular position
theta and beta hat.
The code offers option of comparing solvers in terms of r, theta and total energy calculated, with or
without mass loss, as well as comparing a numerical solver with perturbative expressions for r, theta
or beta hat. For beta hat also an explicit formula can be used for comparison.

## Features
- Calculate orbital parameters, with or without massloss, with Leapfrog integrator or scipy's non-stiff
solvers.
- Calculate orbital parameters beta, x and y numerically, and compare with perturbed solutions
- Comparison plots of energies and orbital parameters

## Requirements
- Python 3.10+
- matplotlib
- numpy
- scipy
- sys

## Installation
1. From zip file or clone repository:

    ```bash 
    git clone https://github.com/wonderbaumer/Project-paper
    ```

2. Navigate into the project folder:
    ```bash
    cd Project-paper
    ```

3. Install requirements:
    ```bash
    pip install -r requirements.txt
    ```

4. Run simulations using main.py
    ```bash
    python main.py
    ```
## Usage
This specifies how to run the main file.
The main file will calculate orbital parameters numerically, using Leapfrog or RK45 and make comparison plots. The comparison plot can be between RK45 and Leapfrog solver, will then be in terms of energy, r or theta, with or without massloss. The comparisons can be between RK45 and perturbed expressions, for r, theta or betahat. User has to either specify a filename for the numerical solutions file they want to use, or run the particle class by specifying time and solver. If solver still not specified, it will revert to default for particle_class. 
Comp_type also has to be specified by user, it can be "r", "theta", "eps_beta", "betahat" or "energy" where each will call different plotting functions from plot.py. 
Note: only some plotting functions are compatible with comparing two numerical solvers.
Note: currently up to first order in epsilon are plotted in default for perturbed expressions, if you want to compare numerical sol without massloss to perturbed expression, set the realistic case of time 1 orbit, t1, and plot only zeroth order of perturbed expressions.


