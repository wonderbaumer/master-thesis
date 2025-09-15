import numpy as np
from forces import gravity
from polar_to_cart import polar_to_cartesian
from leapfrog import leapfrog_algorithm

"""function that calculates the total energy of the particle in orbit,
uses cartesian coordinates"""
def tot_energy(x , y , vx , vy , m , grav):
    """input: x (float), x position
              y (float), y position
              vx (float), x velocity
              vy (float), y velocity 
              grav (float), total gravitational force over x and y dir
              
        returns: tot_E (float), total energy 
        """
    v = np.sqrt(vx**2 + vy**2) #speed
    
    kinetic_energy = 1 / 2 * m * v**2 #formula for kinetic E
    
    pot_energy = grav #formula for potential E
    
    tot_E = kinetic_energy + pot_energy #summing over energies
    
    return tot_E