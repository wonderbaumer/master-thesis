# Project-paper
Particle orbital simulation using numerical solvers and perturbation expressions, all equations are scaled.

## Purpose
This project simulates motion of one particle under the influence of gravitational pull from the Sun,
pressure radiation force (radial component) and mass loss due to sputtering. Simulations are done
with RK4(5) using scipy.integrate.solve_ivp and Leapfrog as 
numerical solvers, and by perturbed expressions for radial distance $\hat{r}$, angular position
$\hat{\theta}$ and $\hat{\beta}$. 
The code offers the option of comparing solvers in terms of $\hat{r}$, $\hat{theta}$ and total energy calculated, with or without mass loss, as well as comparing a numerical solver with perturbative expressions for $\hat{r}$, 
$\hat{\theta}$ or $\hat{\beta}$. For $\hat{\beta}$ also an explicit formula can be used for comparison.
Calculations and plots can also be made for radial and angular velocity, $\hat{v}$ and $\hat{\omega}$.

## Features
- Calculate orbital parameters, with or without massloss, using Leapfrog integrator or scipy's non-stiff
solvers.
- Calculate orbital parameters $\hat{\beta}$, $\hat{r}$ and $\hat{\theta}$ numerically, and compare with perturbed solutions
- Comparison plots of energies and orbital parameters

## Requirements
- Python 3.10+
- matplotlib
- numpy
- scipy

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
This specifies how to run the main file. Running it will make a variety of comparison plots based on input. Not all plots have same attributes. This will explain in that in detail. Note that mandatory input arguments for main file are comp_type and time, and plotting a specified solver will not change labels of the plots.
* "eps_beta": only need mandatory input args to run and provides $\epsilon$ and $\beta$ values corresponding to mass and size range in config file.
* "thetahat": input can be RK4(5) and Leapfrog file, one of the files or a specified solver. If two files input it will plot the comparison of $\hat{\theta}$ between the files, if file or solver it will plot that against the perturbed expression. No combinations will compare a file to solver solution.
* "betahat": input RK4(5) or Leapfrog file, or solver input, rel_fw_err can be True or False. Will compare specified file or solver to perturbed and analytical $\hat{\beta}$ if rel_fw_err==False, else it will plot relative forward error between file or solver and perturbed and analytical $\hat{\beta}$.
* "energy": input RK4(5) and Leapfrog file or solver input, massloss True or False, rel_fw_err True or False. Compares input energies with or without massloss if rel_fw_err ==False, else compares relative forward error in energy between RK4(5) and Leapfrog or solver and file.
* "vhat": input RK4(5), Leapfrog or solver, will compare with perturbed $\hat{v}$.
* "omegahat": input RK4(5), Leapfrog or solver, will compare with perturbed $\hat{\omega}$.


## Bugs
* Changing plot input to specified solver in main does not change labels in plots.
* If choosing file corresponding to mass loss and specifying massloss=False in running main file (or vice versa) where mass loss argument is relevant, the function does not raise error but will run anyways.
* Can not run two specified solvers, need at least one file where applicable. 

## Code Sources and Acknowledgements
Progress bar added to scipy.integrate.solve_ivp
* Lima, T. (2020). Progress Bar with scipy.ipynb [Code]. Github.
https://gist.github.com/thomaslima/d8e795c908f334931354da95acb97e54

Sputtering yields
* Baumann, C., Myrvang, M., & Mann, I. (2020). Dust sputtering within the inner heliosphere: A modelling study [Appendix]. *Annales Geophysicae, 38*, 919-930.

Elemental abundances in the solar wind 
* Killen, R. M., Hurley, D. M., & Farrell, W. M. (2012). The effect on the lunar exosphere of a coronal mass ejection passage [Table] (Vol. 117). *Journal of Geophysical Research*. https://doi.org/10.1029/2011JE004011