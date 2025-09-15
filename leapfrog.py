import numpy as np
from polar_to_cart import polar_to_cartesian
from forces import gravity

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
    
    t = 0 #initial time
    
    x , y , vx , vy = initial_vals[0] #unpack values from nested array
    
    leapfroged_values = [] #initializing list for the calculated values
    
    ax , ay = acc_func(x , y) #unpacking acceleration x and y vals 
    
    while t < t_tot:
        vx_half = vx + dt / 2 * ax #half step calcs for vx
        vy_half = vy + dt / 2 * ay #half step calcs for vy
    
        x = x + dt * vx_half #updating pos x
        y = y + dt * vy_half #updating pos y
        
        ax , ay = acc_func(x , y) #updating acc vals
        
        vx = vx_half + dt / 2 * ax #updating vx
        vy = vy_half + dt / 2 * ay #updating vy
        
        leapfroged_values.append([x , y , vx , vy , ax , ay]) #adding vals to list
        
        t += dt #update time
        
    leapfroged_values = np.array(leapfroged_values)  #list to array
    return leapfroged_values
    

    
    
    


        
    
    