import matplotlib.pyplot as plt
import numpy as np
from forces import *
from polar_to_cart import polar_to_cartesian
from leapfrog import leapfrog_algorithm
from energy import tot_energy

"""plotting the x and y position of the particle"""
def pos_plot(x , y):
    """input: x (float), x position
              y (float), y position
              
        returns: none   """
        
    plt.plot(x , y , label="x and y pos")
    plt.xlabel("x distance (m)")
    plt.ylabel("y distance (m)")
    plt.title("Particle position as function of time")
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
    v0theta = 29.78e3 #initial angular vel in m/s
    
    init_polar = np.array([r0 , theta0 , v0r , v0theta]) #initial values array
    init_cartesian = polar_to_cartesian(init_polar) #initial values to cartesian
    
    dt = 3.16e3 #timestep in s
    t_tot = 3.16e9 #total time in s
    
    pos_and_vel = leapfrog_algorithm(init_cartesian , tot_acc , dt , t_tot) #leapfroging using initial cond
    
    x_pos = pos_and_vel[: , 0] #x pos from leapfrog 
    y_pos = pos_and_vel[: , 1] #y pos from leapfrog
   
    vx = pos_and_vel[: , 2] #vx pos from leapfrog
    vy = pos_and_vel[: , 3] #vy pos from leapfrog
    
    v = np.sqrt(vx**2 + vy**2) #speed 
    
    m = 5.236 * 10**(-7) #mass particle, kg
    
    energy = tot_energy(x_pos , y_pos , vx , vy) #finding energies
    
    energy_plot(t_tot , energy) #plotting energies
    
    #pos_plot(x_pos , y_pos)
    
    
    
    
    
    