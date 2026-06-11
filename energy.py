import numpy as np

"""Calculates the scaled potential and kinetic energies of a particle in orbit, polar coordinates"""
def tot_energy(x , y , vx , vy ,  m , beta , particle_obj):
    """input: x (float), x position
              y (float), y position
              vx (float), x velocity
              vy (float), y velocity
              m (float), mass of particle
              beta (float), beta value
              
        returns: kinetic_energy, pot_energy (tuple), the scaled energies of the particle
        """

    theta_num = np.atan2(y , x) #thetahat
    theta_num = np.unwrap(theta_num) #avoiding discontinuities

    r = np.sqrt(x**2 + y**2)  #radial position
    vr = x * np.cos(theta_num) + y * np.sin(theta_num) #radial velocity
    vtheta = (x * np.sin(theta_num) - y * np.cos(theta_num)) #velocity in theta dir
    
    kinetic_energy = 1 / 2 * m * (vr**2 + vtheta**2) 
    
    pot_energy = -(1 - particle_obj.B * beta) * m / ((1 - particle_obj.B) * r) #scaled potential energy
    
    return kinetic_energy , pot_energy
