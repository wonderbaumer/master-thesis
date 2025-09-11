import numpy as np
from scipy.constants import *

"""scipy constants used:
    G: gravitational constant: 6.6743e-11 m^3kg^-1s^-2"""

m_s = 1.98847e30  #mass of sun, in kg

"""function that converts polar to cartesian coordinates based on 
   arrays containing necessary polar coordinates"""
def polar_to_cartesian(polar_coord):
    """input: polar_coord (array), array containing values for 
              radial position, theta angular position, radial velocity 
              and angular velocity
              
       returns: cartesian_vals (array), array containing cartesian 
                values position x, position y, velocity x and velocity 
                y in same order as they are given in 
    """
    
    cartesian_vals = []  #list for adding cartesian values
    
    for i in polar_coord:
        r , theta , vr , vtheta = polar_coord #separating polar_coord vals
    
        x = r * np.cos(theta) #polar to x position conversion
        y = r * np.sin(theta) #polar to y position conversion
    
        vx = vr * np.cos(theta) - vtheta * np.sin(theta) #conv polar to x vel
        vy = vr * np.sin(theta) + vtheta * np.cos(theta) #conv polar to y vel
    
        cartesian_vals.append(x , y , vx , vy) #adding cartesian coords to list
    
    return np.array(cartesian_vals) #returning array of cartesian list

"""calculates acceleration of the particle in x and y direction
    based on gravitational force between the particle and the sun
    only"""
def gravity(x , y):
    """input: x (float), position in x direction
    
              y (float), position in y direction
              
       returns: acc_x, acc_y, acceleration in x and y direction"""
       
    r = np.sqrt(x**2 + y**2) #radial distance, in km
    acc_x = - G * m_s * x / r**3 #acceleration in x direction,  in kms^-2
    acc_y = - G * m_s * y / r**3 #acceleration in y direction, in kms^-2
    
    return acc_x , acc_y 

"""simple leapfrog algorithm function that uses initial values and
    acceleration from considered forces to calculate position, velocity
    and acceleration of the particle at any given time, in x and y 
    direction"""
def leapfrog_algorithm(initial_vals , acc_func , dt , t_tot):
    """input: initial_vals (array), array containing initial values in
              position x and y and velocity x and y direction
              
              acc_func (function), function calculating acceleration in
              cartesian coordinates based on considered forces, x and y
              direction
              
              dt (float), timestep value
              
              t_tot (float), total simulation time
        
        returns: leapfroged_values (array), array containing position 
        and velocity values, in x and y direction
        """
    
    t = 0
    
    x , y , vx , vy = initial_vals[0:4]
    
    leapfroged_values = []
    
    ax , ay = acc_func(x , y)
    
    while t < t_tot:
        vx_half = vx + dt / 2 * ax
        vy_half = vy + dt / 2 * ay
    
        x = x + dt * vx_half
        y = y + dt * vy_half
        
        ax , ay = acc_func(x , y)
        
        vx = vx_half + dt / 2 * ax
        vy = vy_half + dt / 2 * ay
        
        leapfroged_values.append(x , y , vx , vy)
        leapfroged_values = np.array(leapfroged_values)
        
        t += dt
        
        return leapfroged_values
        


        
    
    