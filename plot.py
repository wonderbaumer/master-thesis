import matplotlib.pyplot as plt
import numpy as np
from forces import *
from polar_to_cart import polar_to_cartesian
from leapfrog import *
from energy import tot_energy
from scipy_solver import *
from constants import *
from analytical_functions import *

"""plotting the x and y position of the particle, two solvers, two subplots"""
def pos_subplot(x_num1 , y_num1 , x_num2 = None , y_num2 = None):
    """input: x_num1 (float), x position from one numerical solver
              y_num1 (float), y position from one numerical solver
              x_num2 (float), x position from another numerical solver
              y_num2 (float), y position from another numerical solver
              
        returns: none   """
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.plot(x_num2 / au , y_num2 / au)
    ax2.plot(x_num1 / au , y_num1 / au)

    ax1.set_aspect("equal")
    ax1.set_xlabel("x distance (AU)")
    ax1.set_ylabel("y distance (AU)")
    ax1.set_title("Runge Kutta 4, t=1000")

    ax2.set_aspect("equal")
    ax2.set_xlabel("x distance (AU)")
    ax2.set_ylabel("y distance (AU)")
    ax2.set_title("Leapfrog, t=1000")
        
    plt.show()

"""plotting numerically calculated beta hat"""
def b_plot(b_num , b_analytical , t):
    """input: betahat (float), beta hat calculated numerically
              t (float), time beta hat has been calculated from
              
       returns: none"""
    fig , ax = plt.subplots()
       
    #time = np.linspace(0 , t / yr , len(b_num)) 
    ax.plot(t , b_num , label = "numerical")
    ax.plot(t , b_analytical , label = "analytical")
    
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("beta hat(t)")
    ax.set_title("Beta hat")
    ax.legend()
    
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    plt.show()

def energy_plot(dt , t_tot , energy1 , energy2):
    kinetic_e1 , _ = energy1
    _ , pot_e1 = energy1
    e_tot1 = kinetic_e1 + pot_e1

    kinetic_e2 , pot_e2 = energy2
    e_tot2 = kinetic_e2 + pot_e2
    
    time = np.linspace(0 , t_tot , 5000) / yr

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    
    ax1.plot(time , kinetic_e1 , label = "Kinetic energy")
    ax1.plot(time , pot_e1 , label = "Potential energy")
    ax1.plot(time , e_tot1 , label = "Total energy")
    
    ax2.plot(time , e_tot2 , label = "Total energy")
    ax2.plot(time , kinetic_e2 , label = "Kinetic energy")
    ax2.plot(time , pot_e2 , label = "Potential energy")
    
    
    #ax1.set_aspect("equal")
    ax1.set_xlabel("Time (yr)")
    ax1.set_ylabel("Energy (J)")
    ax1.set_title("Runge Kutta 4, t=5")
    ax1.legend()

    #ax2.set_aspect("equal")
    ax2.set_xlabel("Time (yr)")
    ax2.set_ylabel("Energy (J)")
    ax2.set_title("Leapfrog, t=5")
    ax2.legend()

    plt.show()

if __name__ == "__main__":
    frog = np.load("leapfrog_numerical.npz")
    x_l = frog["x"][:5000]
    y_l = frog["y"][:5000]
    vx_l = frog["vx"][:5000]
    vy_l = frog["vy"][:5000]

    E_l = tot_energy(x_l , y_l , vx_l , vx_l)
    
    rk45 = np.load("rk45_numerical.npz")
    x_r = rk45["x"][:5000]
    y_r = rk45["y"][:5000]
    vx_r = rk45["vx"][:5000]
    vy_r = rk45["vy"][:5000]

    E_r = tot_energy(x_r , y_r , vx_r , vy_r)

    dt = 3.16e5
    t_tot = 5*3.16e10

    energy_plot(dt , t_tot , E_r , E_l)
    

    #pos_subplot(x_l , y_l , x_num2 = x_r , y_num2 = y_r)
    
    
    