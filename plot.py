import matplotlib.pyplot as plt
import numpy as np
from forces import gravity
from polar_to_cart import polar_to_cartesian
from leapfrog import leapfrog_algorithm
from energy import tot_energy

"""plotting the x and y position of the particle"""
def pos_plot(x , y):
    plt.plot(x , y , label="x and y pos")
    plt.title("Particle position as function of time")
    plt.legend()
    plt.axis("equal")
    plt.show()
    
if __name__ == "__main__":
    
    r0 = 149597871e3 #initial radial dist in m, equals 1 AU
    theta0 = 0 #initial angle in rad, initial position along horizontal
    v0r = 0 #initial radial vel in m/s
    v0theta = 29.78e3 #initial angular vel in m/s
    
    init_polar = np.array([r0 , theta0 , v0r , v0theta]) #initial values array
    init_cartesian = polar_to_cartesian(init_polar) #initial values to cartesian
    
    dt = 3.16e4 #timestep in s
    t_tot = 3.16e8 #total time in s
    
    pos_and_vel = leapfrog_algorithm(init_cartesian , gravity , dt , t_tot) #leapfroging using initial cond
    
    x_pos = pos_and_vel[: , 0]
    y_pos = pos_and_vel[: , 1]
   
    pos_plot(x_pos , y_pos)
    
    
    
    