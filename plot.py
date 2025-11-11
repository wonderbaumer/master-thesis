import matplotlib.pyplot as plt
import numpy as np
from forces import *
from polar_to_cart import polar_to_cartesian
from leapfrog import *
from energy import tot_energy
from scipy_solver import *
from constants import *
from analytical_functions import *

"""plotting the x and y position of the particle"""
def pos_plot(x_num1 , y_num1 , x_num2 = None , y_num2 = None , x_analytical = None , y_analytical = None):
    """input: x_num (float), x position from only numerical solver
              y_num (float), y position from only numerical solver
              x_analytical (float), x pos from analytical expression
              y_analytical (float), y pos from analytical expression
              
        returns: none   """
    
    fig , ax = plt.subplots()
    
    if x_analytical is not None and y_analytical is not None:
        ax.plot(x_analytical / au , y_analytical / au , label = "analytical orbit")
    
    ax.plot(x_num / au , y_num / au , label="orbit")

    
    
    plt.axis("equal")
    plt.xlabel("x distance (AU)")
    plt.ylabel("y distance (AU)")
    plt.title("Particle position as function of time")
    plt.legend()

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

def energy_plot(t , energy):
    kinetic_e , _ = energy
    _ , pot_e = energy
    e_tot = kinetic_e + pot_e
    
    time = np.linspace(0 , t_tot , len(kinetic_e))
    
    plt.plot(time , kinetic_e , label = "kinetic energy")
    plt.plot(time , pot_e , label = "potential energy")
    plt.plot(time , e_tot , label = "total energy")
    
    plt.xlabel("time (s)")
    plt.ylabel("energy (J)")
    plt.title("Energy as function of time")
    plt.legend()
    plt.show()

if __name__ == "__main__":
  
    
    
    