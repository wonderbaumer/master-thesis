import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator as pchip
from config import dat_to_arr, sil_beta , car_beta , size_to_mass
from scipy.constants import c

def inter_func(bval_file):
        size , betaval , _ = dat_to_arr(bval_file)
        size = np.asarray(size).copy()
        size *= 1e-6  ###
        beta_val = pchip(size , betaval , extrapolate = False)

        return beta_val

"""calculates acceleration of the particle in x and y direction
    based on scaled gravitational force between particle and Sun"""
def gravity(x , y , particle_obj):
    """input: x (float), scaled x position
              y (float), scaled y position
              
       returns: acc_x, acc_y (float), scaled acceleration in x and y direction"""
       
    r = np.sqrt(x**2 + y**2) #radial distance
    acc_x = - x / ((1 - particle_obj.B) * r**3) #scaled acceleration in x direction
    acc_y = - y / ((1 - particle_obj.B) * r**3) #scaledacceleration in y direction

    return acc_x , acc_y 

"""calculates mass change from sputtering based on scaled parameters"""
def sputtering(m , epsilon , x , y):
    """input: m (float), scaled mass
        
       return: dmdt (float), mass change as function of time"""

    #dmdt = 0.0
    r = np.sqrt(x**2 + y**2)

    dmdt = - epsilon * m**(2 / 3) * r**(-2)

    return dmdt

"""function that calculates betahat, based on scaled equations"""
def betahat(m , particle_obj):
    """input: m (float), scaled mass of particle

       returns: betahat(float), scaled betahat """
    # b = m**(-1 / 3)
    size = m**(1/3)

    
    r_physical = size * particle_obj.r
    b = particle_obj.beta_func(r_physical) / particle_obj.B

    # r_physmin = 1e-9
    # if (r_physical < r_physmin).any():
    #     raise ValueError(f"Value outside interpolation range")

    return b


"""function that calculates the radial component of the pressure radiation force, 
based on scaled equations"""
def pressure_radial(x , y , m , particle_obj):
    """input: x (float), scaled x position
              y (float), scaled y position
              m (float), scaled mass of particle

        returns: ax , ay (tuple), acceleration in x and y dir for pressure radiation force"""
    
    r = np.sqrt(x**2 + y**2) #scaled radial distance

    ax = x * betahat(m , particle_obj) * particle_obj.B / ((1 - particle_obj.B) * r**3) #scaled acceleration in x dir
    ay = y * betahat(m , particle_obj) * particle_obj.B / ((1 - particle_obj.B) * r**3) #scaled acceleration in y dir

    return ax , ay

"""calculates Poynting-Robertson drag based on scaled parameters"""
def pr_drag(x , y , vx , vy , m , particle_obj):
    """input: x (float), scaled x position
              y (float), scaled y position
              vx (float), scaled x velocity
              vy (float), scaled y velocity
              m (float), scaled particle mass

        returns ax , ay (tuple), acceleration in x and y direction"""
    
    r = np.sqrt(x**2 + y**2) #scaled radial distance

    theta = np.atan2(y , x)

    A = -betahat(m , particle_obj) * particle_obj.B * particle_obj.delta / ((1 - particle_obj.B) * r**3)
    #A = -betahat(m , particle_obj) * particle_obj.B * particle_obj.delta / ((1 - particle_obj.B) * r**3)
    x_dir = 2 * np.cos(theta) * (x * vx + y * vy) - np.sin(theta) * (x * vy - y * vx)
    y_dir = 2 * np.sin(theta) * (x * vx + y * vy) + np.cos(theta) * (x * vy - y * vx)

    ax = A * x_dir
    ay = A * y_dir

    return ax , ay

"""function that calculates total acceleration based on pressure radiation force, gravity and Poynting-Robertson drag"""
def tot_acc(x , y , vx , vy , m , particle_obj):
    """input: x (float), scaled x position
              y (float), scaled y position
              vx (float), scaled x velocity
              vy (float), scaled y velocity
              m (float), scaled mass of particle

        returns: ax , ay (array), scaled acceleration in x and y dir"""
    
    px , py = pressure_radial(x , y , m , particle_obj) #decomposing pressure radiation force
    gx , gy = gravity(x , y , particle_obj) #decomposing gravitational force
    prx , pry = pr_drag(x , y , vx , vy , m , particle_obj) #decomposing Poynting-Robertson force

    ax = px + gx + prx  #total acceleration in x dir
    ay = py + gy + pry  #total acceleration in y dir 
    
    return ax , ay

    
    