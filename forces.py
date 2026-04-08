from scipy.constants import G , c
import numpy as np
from config import S_s , q_pr , m_s , rho_s , rho_c , au

"""calculates acceleration of the particle in x and y direction
    based on gravitational force between the particle and the sun
    only"""
def gravity(x , y):
    """input: x (float), position in x direction in m
    
              y (float), position in y direction in m
              
       returns: acc_x, acc_y (float), acceleration in x and y direction"""
       
    r = np.sqrt(x**2 + y**2) #radial distance, in m
    acc_x = - G * m_s * x / r**3 #acceleration in x direction,  in ms^-2
    acc_y = - G * m_s * y / r**3 #acceleration in y direction, in ms^-2
    
    return acc_x , acc_y 

"""calculates mass change from sputtering"""
def sputtering(m , epsilon):
    """input: m (float), mass in kg
              epsilon (float), epsilon value 
        
       return: dmdt (float), mass change as function of time"""
    
    dmdt = - epsilon * m**(2 / 3) #total mass change
    
    return dmdt

"""calculates radius of particle based on mass"""
def radius(m):
    """input: m (float), mass of particle in kg
       
       returns: r (float), radius of particle in m"""
       
    r = (3 * m / (4 * rho_s * np.pi))**(1 / 3) #radius based on perfect sphere assumptions, in m SILICATE ONLY
    
    return r

"""function that calculates the radial component of the radiation pressure
force from the solar wind hitting the particle"""
def pressure_radial(x , y , m):
    """input: x (float), cartesian x coordinate for position in m
              y (float), cartesian y coordinate for pos in m
              m (float), mass of particle, in kg

        returns: pressure_force_rad (float), radiation pressure force in N, rad comp"""
        
    r = np.sqrt(x**2 + y**2) #radial distance of particle from Sun in m
    r_par = radius(m) #radius of the particle in m
    
    A = np.pi * r_par**2 #cross section area of particle in m^2
    
    s = S_s * (au / r)**2 #radiation flux density at distance r from Sun, in Wm^-2
    
    pressure_force_rad = s * A * q_pr / c #formula pressure radiation force, in N
    
    return pressure_force_rad

"""function that calculates beta, the ratio between pressure radiation force
and gravity"""
def beta(x , y , m):
    """input: x (float), cartesian x coordinate for position in m
              y (float), cartesian y coordinate for pos in m
              m (float), mass of particle in kg
        
       returns: b (float), ratio between radiation and gravitation force"""
    
    gx , gy = gravity(x , y) #gravitational acceleration in x and y dir in ms^-2
    g_abs = np.sqrt(gx**2 + gy**2) #absolute value gravity
    
    Frad = pressure_radial(x , y , m) #in N
    arad = Frad / m #radiation acc in ms^-2
    
    b = arad / g_abs #ratio radiation pressure acc to gravity
    
    return b

"""calculates acceleration due to Poynting-Robertson drag"""
def pr_drag(x , y , vx , vy , m):
    """input: x (float), cartesian x coordinate for position in m
              y (float), cartesian y coordinate for pos in m
              vx (float), cartesian x velocity in ms^-1
              vy (float), cartesian y velocity in ms^-1
              m (float), mass of particle, in kg

        returns: ax , ay (tuple), acceleration in x and y direction"""
        
    r = np.sqrt(x**2 + y**2) #radial distance of particle from Sun in m
    theta = np.atan2(y , x)

    b = beta(x , y , m)

    grav_term = G * m_s / r**4

    xvel_terms = -2 * vx * (x**2 + y**2) - x * y * vy #x velocity terms
    yvel_terms = -x * y * vx - 2 * vy * (y**2 + x**2) #y velocity terms

    ax = grav_term * b / c * xvel_terms #acceleration in x direction
    ay = grav_term * b / c * yvel_terms #acceleration in y direction

    return ax , ay

"""function that calculates total acceleration given radiation pressure force
and gravitational force only"""
def tot_acc(x , y , vx , vy , m):
    """input: x (float), cartesian x coordinate for position in m
              y (float), cartesian y coordinate for pos in m
              m (float), mass of particle in kg
              
        returns: ax , ay (array), calculated acceleration in x and y dir, in ms^-2"""
        
    pressure_acc = pressure_radial(x , y , m) / m #acceleration due to pressure radiation, in ms^-2
    gx , gy = gravity(x , y) #gravitational acceleration in x and y dir, in ms^-2
    
    b = pressure_acc / np.sqrt(gx**2 + gy**2) #beta calcs
    prx , pry = pr_drag(x , y , vx , vy , m) 
    
    #ax = gx * (1 - b) #acc in x dir in ms^-2 without drag
    #ay = gy * (1 - b) #acc in y dir in ms^-2 without drag

    ax = gx * (1 - b) + prx #acc in x dir in ms^-2 with drag
    ay = gy * (1 - b) + pry #acc in y dir in ms^-2 with drag
    
    return ax , ay



