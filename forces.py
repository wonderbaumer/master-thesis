from scipy.constants import *
import numpy as np
from polar_to_cart import polar_to_cartesian
from cst_table import *
from constants import *

"""scipy constants used:
    G: gravitational constant: 6.6743e-11 m^3kg^-1s^-2
    c: speed of light in vacuum: 299792458 ms^-1
    
    r_par, m_par calculated on paper, necessary for calcs in this py, but not
    in class py"""

r_par = 500e-9
m_par = 1.30899694e-15 #mass particle, kg

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
def pressure_radial(x , y , r_par):
    """input: x (float), cartesian x coordinate for position
              y (float), cartesian y coordinate for pos
              r_par (float), radius of particle
        returns: pressure_force_rad (float), radiation pressure force, rad comp"""
        
    r = np.sqrt(x**2 + y**2) #radial distance of particle from Sun
    A = np.pi * r_par**2 #cross section area of particle
    
    s = S_s * (au / r)**2 #radiation flux density at distance r from Sun
    
    pressure_force_rad = s * A * q_pr / c #formula pressure radiation force
    
    return pressure_force_rad

"""function that calculates beta, ratio between pressure radiation force
and gravity"""
def beta(x , y , r_par):
    """input: x (float), cartesian x coordinate for position
              y (float), cartesian y coordinate for pos
              r_par (float), radius of particle
        
       returns: b (float), ratio between rad and g force"""
    
    gx , gy = gravity(x , y) #gravitational acceleration in x and y dir
    g_abs = np.sqrt(gx**2 + gy**2) #absolute value gravity
    
    Frad = pressure_radial(x , y , r_par)
    arad = Frad / m_par #radiation acc
    
    b = arad / g_abs #ratio radiation pressure acc to gravity
    
    return b

"""function that calculates total acceleration given radiation pressure force
and gravitational force only"""
def tot_acc(x , y , beta):
    """input: x (float), cartesian x coordinate for position
              y (float), cartesian y coordinate for pos
              beta (float), Fr/Fg
              
        returns: ax , ay (array), lists of new acceleration in x and y dir"""
        
    gx , gy = gravity(x , y) #gravitational acceleration in x and y dir
    
    ax = gx * (1 - beta) #acc in x dir
    ay = gy * (1 - beta) #acc in y dir
    
    return ax , ay
   

if __name__ == "__main__":
    init_polar = np.array([r0 , theta0 , v0r , v0theta]) #initial values array
    init_cartesian = polar_to_cartesian(init_polar) #initial values to cartesian
    x = init_cartesian[0][0]
    y = init_cartesian[0][1]
    
    
    
    