# Project-paper
Particle orbital simulation using numerical solvers and perturbation expressions, all equations are scaled.

## Purpose
This project simulates motion of one particle under the influence of gravitational pull from the Sun,
 radiation pressure force, Poynting Robertson drag from solar photons, and mass loss due to sputtering. Furthermore, the complete set of values for the radiation pressure-to-gravity ratio, $\beta$, has been provided, denoted numerical $\beta$. Simulations are done with RK4(5) using scipy.integrate.solve_ivp as numerical solver, both for the reference solution as well as the perturbed solutions. The solver validatation is done by considering only the gravitational pull and radiation pressure force, and evaluating energy conservation. The perturbation has been done using a fast timescale (changes over one orbit) and slow timescale (changes over several orbits), and the variables have been perturbed in $\epsilon$, which is the mass loss rate due to sputtering. The Poynting-Robertson contribution to perturbation is handled by writing $\delta$, the drag term, as function of $\epsilon$. An analytical approximation for $\beta$ is used, denoted analytical $\beta$. The perturbed solutions are validated by comparing with the numerical solutions, where both use the analytical $\beta$.

This code has the options of comparing numerical solvers in terms of energy considering conservative effects and comparing numerical to perturbed solutions for $\hat{r}$, $\hat{\omega}$ and $\hat{\beta}$. The numerical solver can evaluate the numerical $\beta$ values and compare with the results using the analytical $\beta$ values, which allow the user to investigate the overestimation the analytical $\beta$ does on the orbit. Furthermore, the calculations can be done both for silicate (MgFeSiO_4) and carbon (C), for various sizes. Using this feature one can also investigate the differences in orbit between silicate and carbon under the same initial conditions.

## Features
- Calculate orbital parameters using scipy's non-stiff solver RK4(5).
- Calculate orbital parameters $\hat{\beta}$, $\hat{r}$ and $\hat{\theta}$ numerically, and compare with perturbed solutions.
- Comparison plots of energies and orbital parameters
- Complete dictionary of sputtering yields for carbon and silicate, for nine solar wind ions, and for fast and slow solar wind, as well as coronal mass ejection (CME) conditions.
- Dictionaries of perturbed and numerical lifetimes which can be used as point of comparison of lifetimes between numerical, perturbed solutions, and theoretical values. 
- $\beta$ calculated from Mie theory

## Requirements
- Python 3.10+
- matplotlib
- numpy
- scipy
- sympy
- tqdm

## Installation
1. From zip file or clone repository:

    ```bash 
    git clone https://github.com/wonderbaumer/master-thesis
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
All plots produced for the master thesis are produced through the main file. The main function initiates plotting files depending on user input. This is run for all the plots in the file and is hashed out initially, so depending on which plot the user want to reproduce one has to un-hash the selected lines. Specified for each plot, if applicable, is the file necessary to make the plot. Those are made through the make_file function in the main file, and all necessary files are made and hashed out initially, so to make the desired files one has to un-hash selected lines.

If one want to run simulations outside the main file, the numerical simulations are performed through particle_class file and perturbed simulations in pert_variable_eps file. Into the functions one need to pass arguments of dust properties, which are defined in dust_properties. The simulation results can be plotted in plot file, but please note the plot functions expect a .npz file for the numerical plots.


## Bugs


## Code Sources and Acknowledgements
Progress bar added to scipy.integrate.solve_ivp
* Lima, T. (2020). Progress Bar with scipy.ipynb [Code]. Github.
https://gist.github.com/thomaslima/d8e795c908f334931354da95acb97e54

Sputtering yields
* Baumann, C., Myrvang, M., & Mann, I. (2020). Dust sputtering within the inner heliosphere: A modelling study [Appendix]. *Annales Geophysicae, 38*, 919-930.

Elemental abundances in the solar wind 
* Killen, R. M., Hurley, D. M., & Farrell, W. M. (2012). The effect on the lunar exosphere of a coronal mass ejection passage [Table] (Vol. 117). *Journal of Geophysical Research*. https://doi.org/10.1029/2011JE004011

Arrow in epsilon and delta vs B plot
* User Dietrich answer to question "Draw arrow outside plot in Matplotlib" (2014) [Code]. https://stackoverflow.com/questions/23922804/draw-arrow-outside-plot-in-matplotlib

Calculated $\beta$ values in files ac_radpr_prdrag_sun1au.dat and sil_radpr_prdrag_sun1au.dat
* Li, A. (personal communication, March 1, 2026) [Table].
