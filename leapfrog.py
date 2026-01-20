import numpy as np 
from config import mhat0 , t5 , t6 , t7 , init_cart_scaled
from forces_scaled import tot_acc, sputtering, betahat

"""simple leapfrog algorithm that usesinitial values and acceleration from considered forces to 
    calculate position, velocity, mass and beta value of the particle at any given time, in x and y 
    direction. All parameters are scaled."""
def leapfrog_algorithm(initial_vals , acc_func , time , massloss = None):
    """input: initial_vals (array), array containing scaled cartesian x, y, vx, vy, m.
              
              acc_func (function), function calculating total scaled acceleration in x and y dir
              
              dt (float), scaled timestep value
              
              t_tot (float), total simulation time in s

              massloss (str), default: None, else method of massloss must be defined 
        
        returns: leapfroged_values (array) , b_vals (array), array containing position, velocity, time 
        and mass, array containing beta values
        """
    
    x , y , vx , vy = initial_vals #unpack values from nested array

    dt , t_tot = time #unpack time values
    t = 0 #initial time

    mhat = mhat0 #initial scaled mass
    bhat = betahat(mhat) #initial scaled beta

    N = int(t_tot / dt) + 1 #number of timesteps
    lf_vals = np.zeros((N, 6)) #array to store leapfrogged values
    lf_vals[0] = [x, y, vx, vy, mhat, bhat]

    ax , ay = acc_func(x , y , mhat) #unpacking acceleration x and y vals 

    vx_half = vx + 0.5 * dt * ax #half-stepping x velocity
    vy_half = vy + 0.5 * dt * ay #half-stepping y velocity

    #half-stepping mass calcs if massloss is considered
    if massloss is not None:
        m_half = mhat + 0.5 * dt * massloss(mhat) 

    for i in range(1, N):
        x += dt * vx_half #pos x calcs
        y += dt * vy_half #pos y calcs

        #updating mass if massloss is considered
        if massloss is not None:
            mhat = m_half + 0.5 * dt * massloss(m_half) 
            m_half += dt * massloss(m_half) #updating mass

        bhat = betahat(mhat)
            
        ax , ay = acc_func(x , y , mhat) #acceleration calcs
        
        vx_half += dt * ax #updating vx_half
        vy_half += dt * ay #updating vy_half

        vx = vx_half - 0.5 * dt * ax #updating vx
        vy = vy_half - 0.5 * dt * ay #updating vy

        lf_vals[i] = [x, y, vx, vy, mhat, bhat] #leapfroged_values 
    
    return lf_vals 
    
if __name__ == "__main__":
    dt , t_tot = t6 
    
    lf_vals  = leapfrog_algorithm(init_cart_scaled , tot_acc , t6)
    x , y , vx , vy , m , beta = lf_vals[: , 0] , lf_vals[: , 1] , lf_vals[: , 2] , lf_vals[: , 3] , lf_vals[: , 4] , lf_vals[: , 5]
    print(np.sqrt(x**2 + y**2))
    #np.savez(f"C:/Users/cecil/Documents/Project-paper/Files/leapfrog_t6_masslossFalse_scaledeqs" , x = x , y = y , vx = vx , vy = vy , m = m , b = b_vals)