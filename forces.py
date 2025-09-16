from scipy.constants import *
import numpy as np

"""scipy constants used:
    G: gravitational constant: 6.6743e-11 m^3kg^-1s^-2
    c: speed of light in vacuum: 299792458 ms^-1"""

m_s = 1.98847e30  #mass of sun, in kg
s_s = 1361 #solar constant, in Wm^-2
q_pr = 1 #radiation pressure coefficient, unitless

"""calculates acceleration of the particle in x and y direction
    based on gravitational force between the particle and the sun
    only"""
def gravity(x , y):
    """input: x (float), position in x direction
    
              y (float), position in y direction
              
       returns: acc_x, acc_y (float), acceleration in x and y direction"""
       
    r = np.sqrt(x**2 + y**2) #radial distance, in m
    acc_x = - G * m_s * x / r**3 #acceleration in x direction,  in ms^-2
    acc_y = - G * m_s * y / r**3 #acceleration in y direction, in ms^-2
    
    return acc_x , acc_y 

"""function that calculates the radial part of the radiation pressure
force from the solar wind hitting the particle"""
def pressure_radial(x , y , r0 , r_par):
    """input: x (float), cartesian x coordinate for position
              y (float), cartesian y coordinate for pos
              r0 (float), initial distance between par and Sun
              r_par (float), radius of particle
        returns: pressure_force_rad (float), radiation pressure force, rad comp"""
        
    r = np.sqrt(x**2 + y**2)
    A = np.pi * r_par**2
    
    s = s_s * (r0 / r)**2
    
    pressure_force_rad = s * A * q_pr / c
    
    return pressure_force_rad


#def acc_beta() 
    

    
    