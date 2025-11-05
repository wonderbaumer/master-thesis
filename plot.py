import matplotlib.pyplot as plt
import numpy as np
from forces import *
from polar_to_cart import polar_to_cartesian
from leapfrog import leapfrog_algorithm
from energy import tot_energy
from scipy_solver import *
from constants import *

"""plotting the x and y position of the particle"""
def pos_plot(x , y):
    """input: x (float), x position
              y (float), y position
              
        returns: none   """
        
    plt.plot(x / au , y / au , label="x and y pos")
    plt.axis("equal")
    plt.xlabel("x distance (AU)")
    plt.ylabel("y distance (AU)")
    plt.title("Particle position as function of time")
    plt.legend()

    plt.show()

"""plotting numerically calculated beta hat"""
def betahat_plot(betahat , t):
    """input: betahat (float), beta hat calculated numerically
              t (float), time beta hat has been calculated from
              
       returns: none"""
       
    time = np.linspace(0 , t / yr , len(betahat)) 
    plt.plot(time , betahat)
    
    plt.xlabel("Time (yr)")
    plt.ylabel("beta hat(t)")
    plt.title("Beta hat numerical")
    plt.legend()
    plt.show()

def energy_plot(t , energy):
    kinetic_e , _ = energy
    _ , pot_e = energy
    e_tot = kinetic_e + pot_e
    
    time = np.linspace(0 , t_tot , len(kinetic_e))
    
    plt.plot(time , kinetic_e , label = "kinetic energy")
    plt.plot(time , pot_e , label = "potential energy")
    plt.plot(time , e_tot , label = "total energy")
    
    plt.xlabel("time (s)")
    plt.ylabel("energy (J)")
    plt.title("Energy as function of time")
    plt.legend()
    plt.show()
    
if __name__ == "__main__":
    
    theta0 = 0 #initial angle in rad, initial position along horizontal
    v0r = 0 #initial radial vel in m/s
    v0theta = 22085 #initial angular vel in m/s
    
    init_polar = np.array([r0 , theta0 , v0r , v0theta]) #initial values array
    init_cartesian = polar_to_cartesian(init_polar) #initial values to cartesian
    
    dt = 3.16e2 #timestep in s
    t0 = 0 #initial time in s
    t_tot = yr #total time in s
    t_span = (t0 , t_tot) #tuple of start and end time

    #pos_scipy = particle_motion(pos_vel , t_span , init_cartesian[0])
    
    betas = leapfrog_algorithm(init_cartesian , tot_acc , sputtering , dt , t_tot) #leapfroging using initial cond
    betahat = betas / beta_0
    
    #betahat_plot(betahat , t_tot)
    
    #x_pos = pos_and_vel[: , 0] #x pos from leapfrog 
    #y_pos = pos_and_vel[: , 1] #y pos from leapfrog
   
    #vx = pos_and_vel[: , 2] #vx pos from leapfrog
    #vy = pos_and_vel[: , 3] #vy pos from leapfrog
    
    #x_scipy = pos_scipy.y[0]
    #y_scipy = pos_scipy.y[1]
    
    #vx_scipy = pos_scipy.y[2]
    #vy_scipy = pos_scipy.y[3]
    
    #v = np.sqrt(vx**2 + vy**2) #speed 
    
    #energy = tot_energy(x_scipy , y_scipy , vx_scipy , vy_scipy) #finding energies
    
    #energy_plot(t_tot , energy) #plotting energies
    
    #pos_plot(x_pos , y_pos)
    
    
    
    
    
    