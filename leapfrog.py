import numpy as np
from forces import * 
from config import *

"""simple leapfrog algorithm that uses initial values and acceleration from considered forces to 
    calculate position, velocity and acceleration of the particle at any given time, in x and y 
    direction"""
def leapfrog_algorithm(initial_vals , acc_func , dt , t_tot , massloss = None):
    """input: initial_vals (array), array containing cartesian initial position and velocity, x and y
              in m.
              
              acc_func (function), function calculating total acceleration in
              cartesian coordinates based on considered forces, x and y
              direction
              
              dt (float), timestep value in s
              
              t_tot (float), total simulation time in s

              massloss (str), default: None, else method of massloss must be defined 
        
        returns: leapfroged_values (array) , b_vals (array), array containing position, velocity, time 
        and mass, array containing beta values
        """
    
    x , y , vx , vy = initial_vals #unpack values from nested array
    
    leapfroged_values = [] #initializing list for the calculated values
    
    
    t = 0 #initial time
    m = m_par #initial mass
    b = beta(x , y , m) #initial beta 
    b_vals = [b] #initializing list of beta values

    ax , ay = acc_func(x , y , m) #unpacking acceleration x and y vals 

    vx_half = vx + 0.5 * dt * ax #half-stepping x velocity
    vy_half = vy + 0.5 * dt * ay #half-stepping y velocity

    #half-stepping mass calcs if massloss is considered
    if massloss is not None:
        m_half = m + 0.5 * dt * massloss(m) 

    while t < t_tot:
        x += dt * vx_half #pos x calcs
        y += dt * vy_half #pos y calcs

        #updating mass if massloss is considered
        if massloss is not None:
            m = m_half + 0.5 * dt * massloss(m_half)  #mass calcs
            m_half += dt * massloss(m_half) #updating mass

        b_vals.append(beta(x, y, m)) #beta calcs
            
        ax , ay = acc_func(x , y , m) #acceleration calcs
        
        vx_half += dt * ax #updating vx_half
        vy_half += dt * ay #updating vy_half

        vx = vx_half - 0.5 * dt * ax #updating vx
        vy = vy_half - 0.5 * dt * ay #updating vy

        leapfroged_values.append([x , y , vx , vy , m]) #putting values into list
        
        t += dt #update time
        
    leapfroged_values = np.array(leapfroged_values)  #leapfrog values list to array
    b_vals = np.array(b_vals[:-1]) #beta values list to array, drop final value for same length
    
    return leapfroged_values , b_vals
    
if __name__ == "__main__":
    dt , t_tot = t4   
    lf_vals , b_vals = leapfrog_algorithm(init_cartesian , tot_acc , dt , t_tot)
    x , y , vx , vy , m = lf_vals[: , 0] , lf_vals[: , 1] , lf_vals[: , 2] , lf_vals[: , 3] , lf_vals[: , 4]