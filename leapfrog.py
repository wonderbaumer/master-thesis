import numpy as np
from polar_to_cart import polar_to_cartesian
from forces import *

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
    
    ax , ay = acc_func(x , y , r_par) #unpacking acceleration x and y vals 
    
    while t < t_tot:
        vx_half = vx + dt / 2 * ax #half step calcs for vx
        vy_half = vy + dt / 2 * ay #half step calcs for vy
    
        x = x + dt * vx_half #updating pos x
        y = y + dt * vy_half #updating pos y
        
        ax , ay = acc_func(x , y , r_par) #updating acc vals
        
        vx = vx_half + dt / 2 * ax #updating vx
        vy = vy_half + dt / 2 * ay #updating vy
        
        leapfroged_values.append([x , y , vx , vy , ax , ay]) #adding vals to list
        
        t += dt #update time
        
    leapfroged_values = np.array(leapfroged_values)  #list to array
    return leapfroged_values
    

if __name__ == "__main__":  
    
    theta0 = 0 #initial angle in rad, initial position along horizontal
    v0r = 0 #initial radial vel in m/s
    v0theta = 29.78e3 #initial angular vel in m/s
    
    init_polar = np.array([r0 , theta0 , v0r , v0theta]) #initial values array
    init_cartesian = polar_to_cartesian(init_polar) #initial values to cartesian
    
    dt = 3.16e5 #timestep in s
    t_tot = 3.16e10 #total time in s
    
    pos_and_vel = leapfrog_algorithm(init_cartesian , tot_acc , dt , t_tot) #leapfroging using initial cond
    x_pos = pos_and_vel[: , 0] #x pos from leapfrog 
    y_pos = pos_and_vel[: , 1] #y pos from leapfrog
    
    print(beta(x_pos, y_pos, r_par))
    
    
    
    
   
    
    
    
    


        
    
    