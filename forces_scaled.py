from scipy.constants import *
import numpy as np
import sys
sys.path.insert(1, 'C:/Users/Cecilie.Bamer/Documents/Project-paper/')
from constants import *
from config import *

"""scipy constants used:
    G: gravitational constant: 6.6743e-11 m^3kg^-1s^-2
    c: speed of light in vacuum: 299792458 ms^-1
    
    r_par, m_par are hand-calculated"""

"""calculates acceleration of the particle in x and y direction
    based on gravitational force between the particle and the sun
    only"""
def gravity(x , y):
    """input: x (float), position in x direction in m
    
              y (float), position in y direction in m
              
       returns: acc_x, acc_y (float), acceleration in x and y direction"""
       
    r = np.sqrt(x**2 + y**2) #radial distance
    acc_x = - x / ((1 - beta0) * r**3) #acceleration in x direction
    acc_y = - y / ((1 - beta0) * r**3) #acceleration in y direction

    return acc_x , acc_y 

def mhat(m):
    mhat = m / m_par

    return mhat

"""calculates mass change from sputtering"""
def sputtering(m):
    """input: m (float), mass in kg
        
       return: dmdt (float), mass change as function of time"""

    #m_hat = mhat(m) 
    
    #dmdt = - eps(m) * m_hat**(2 / 3)

    dmdt = - eps() * m**(2 / 3)
    
    return dmdt

"""function that calculates beta, the ratio between pressure radiation force
and gravity"""
def betahat(m):
    """input: x (float), cartesian x coordinate for position in m
              y (float), cartesian y coordinate for pos in m
              m (float), mass of particle in kg
        
       returns: b (float), ratio between radiation and gravitation force"""
    #m_hat = mhat(m)
    
    #betahat = m_hat**(-1 / 3)
    betahat = m**(-1 / 3)
    
    return betahat

"""function that calculates the radial component of the radiation pressure
force from the solar wind hitting the particle"""
def pressure_radial(x , y , m):
    """input: x (float), cartesian x coordinate for position in m
              y (float), cartesian y coordinate for pos in m
              m (float), mass of particle, in kg

        returns: pressure_force_rad (float), radiation pressure force in N, rad comp"""
    
    r = np.sqrt(x**2 + y**2)

    ax = x * betahat(m) * beta0 / ((1 - beta0) * r**3)
    ay = y * betahat(m) * beta0 / ((1 - beta0) * r**3)

    return ax , ay

"""function that calculates total acceleration given radiation pressure force
and gravitational force only"""
def tot_acc(x , y , m):
    """input: x (float), cartesian x coordinate for position in m
              y (float), cartesian y coordinate for pos in m
              m (float), mass of particle in kg
              
        returns: ax , ay (array), calculated acceleration in x and y dir, in ms^-2"""
    
    px , py = pressure_radial(x , y , m)
    gx , gy = gravity(x , y)

    ax = px + gx
    ay = py + gy    
    
    return ax , ay

if __name__ == "__main__":
    x , y , vx , vy = init_cart_scaled

    print(tot_acc(x , y , m_par / m_par))
    