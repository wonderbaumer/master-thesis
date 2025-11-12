import numpy as np
from polar_to_cart import polar_to_cartesian
from forces import *
from analytical_functions import *
from constants import * 
import matplotlib.pyplot as plt

"""simple leapfrog algorithm function that uses initial values and
    acceleration from considered forces to calculate position, velocity
    and acceleration of the particle at any given time, in x and y 
    direction"""
def leapfrog_algorithm(initial_vals , acc_func , dt , t_tot , massloss = None , scaled = False):
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
    x , y , vx , vy = initial_vals #unpack values from nested array
    
    leapfroged_values = [] #initializing list for the calculated values
    b_vals = []
    
    t = 0 #initial time
    
    m = m_par
    b = beta(x , y , m)
    b_vals.append(b)

    if scaled:
        x , y  = x / X , y / X
        vx , vy =  vx / V , vy / V

        ax , ay = acc_func(x , y , t)
        
    else:
        ax , ay = acc_func(x , y , m) #unpacking acceleration x and y vals 
        
    while t < t_tot:
        vx_half = vx + dt / 2 * ax #half step calcs for vx
        vy_half = vy + dt / 2 * ay #half step calcs for vy
        
        if massloss is not None:
            dmdt = massloss(m)
            m_mid = m + 0.5 * dt * dmdt
            
            """correcting for numerical instabilities in mass cals"""
            dm_mid = massloss(m_mid)
            m += dm_mid * dt
            
            b = beta(x , y , m)
            b_vals.append(b)
           
        x = x + dt * vx_half #updating pos x
        y = y + dt * vy_half #updating pos y
        
        if scaled:
            xhat , yhat = x / X , y / X
            vxhat , vyhat = vx / V , vy / V
            ax , ay = acc_func(x , y , t)
            
        else:
            ax , ay = acc_func(x , y , m)
        
        vx = vx_half + dt / 2 * ax #updating vx
        vy = vy_half + dt / 2 * ay #updating vy
            
        leapfroged_values.append([x , y , vx , vy , ax , ay]) #adding vals to list
        
        t += dt #update time
      
    if scaled:
        leapfroged_values = [(X * xi , X * yi , V * vx , V * vy) for xi , yi ,
                             vx , vy , axi , ayi in leapfroged_values]
        
    leapfroged_values = np.array(leapfroged_values)  #list to array
    b_vals = np.array(b_vals)
    
    
    return leapfroged_values , b_vals
    

if __name__ == "__main__":
    v0theta = 2.19013101e+04 #initial angular vel in m/s

    init_polar = np.array([r0 , theta0 , v0r , v0theta]) #initial values array
    init_cartesian = polar_to_cartesian(init_polar) #initial values to cartesian

    dt = 3.16e5 #timestep in s
    t_tot = 5*3.16e10#total time in s

    scaled_dt = dt / T
    scaled_ttot = t_tot / T

    #orbit_a0 = leapfrog_algorithm(init_cartesian , acceleration , scaled_dt , scaled_ttot 
                                #, scaled = True)
    orbit_n , _ = leapfrog_algorithm(init_cartesian, tot_acc , dt , t_tot , sputtering)
    x = orbit_n[: , 0]
    y = orbit_n[: , 1]
    vx = orbit_n[: , 2]
    vy = orbit_n[: , 3]
    ax = orbit_n[: , 4]
    ay = orbit_n[: , 5]
    print(len(ay))
    print(np.sqrt(x**2+y**2) / au)

  
    
    
    
    


    
    
    
    
    
   
    
    
    
    


        
    
    