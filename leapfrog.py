import numpy as np 
import time
from tqdm import tqdm
from config import t6 , t7 , t8 , t9 , t10
from forces_scaled import tot_acc , sputtering, betahat

"""Simple leapfrog algorithm that use acceleration to calculate position, velocity, 
   mass and beta value of the particle at any given time, in x and y direction. 
   All parameters are scaled."""
def leapfrog_algorithm(initial_vals , acc_func , time , particle_obj , epsilon = None , massloss = None 
                       , drag = True):
    """input: initial_vals (array), array containing scaled cartesian x, y, vx, vy, m.
              acc_func (function), function calculating total scaled acceleration in x and y dir
              time (tuple), consisting of dt, t_tot, timestep and total simulation time
              particle_obj, instance containing simulation information
              epsilon (float), mass loss rate value, default:None 
              massloss (str), default: None, else method of massloss must be defined 
        
        returns: leapfroged_values (array), containing position, velocity, mass, beta, time
        """
    
    x , y , vx , vy , mhat = initial_vals #unpack initial values

    dt , t_tot = time #unpack time values
    
    bhat = 1.0 #initial scaled beta

    N = round(t_tot / dt) + 1 #number of timesteps
    lf_vals = np.zeros((N , 7)) #array to store leapfrogged values
    dt = t_tot / (N - 1)

    t = 0.0

    pbar = tqdm(total = N)

    lf_vals[0] = [x , y , vx , vy , mhat , bhat , t]

    ax , ay = acc_func(x , y , vx , vy , mhat , particle_obj , drag) #unpacking acceleration x and y vals 

    vx_half = vx + 0.5 * dt * ax #half-stepping x velocity
    vy_half = vy + 0.5 * dt * ay #half-stepping y velocity

    for i in range(1 , N):
        t += dt
        x += dt * vx_half #pos x calcs
        y += dt * vy_half #pos y calcs

        #updating mass if massloss is considered
        if massloss is not None:
            mhat += 0.5 * dt * massloss(mhat , epsilon , x , y) 

        bhat = betahat(mhat , particle_obj)
            
        ax , ay = acc_func(x , y , vx_half , vy_half , mhat , particle_obj , drag) #acceleration calcs
        
        vx_half += dt * ax #updating vx_half
        vy_half += dt * ay #updating vy_half

        lf_vals[i] = [x , y , vx_half , vy_half , mhat , bhat , t] #leapfrogged values 

        pbar.update(1)

    pbar.close()
    
    return lf_vals 

    
