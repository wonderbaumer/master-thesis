import numpy as np
from scipy.constants import *
from constants import *
from config import *

"""function that calculates the total energy of the particle in orbit,
uses cartesian coordinates"""
def tot_energy(x , y , vx , vy , m , beta = beta0):
    """input: x (float), x position
              y (float), y position
              vx (float), x velocity
              vy (float), y velocity
              m (float), mass of particle
              beta (float), default: beta0 if massloss not considered, else user-specified
              
        returns: kinetic_energy, pot_energy (tuple), energies
        """
        
    v = np.sqrt(vx**2 + vy**2) #speed
    
    kinetic_energy = 1 / 2 * m * v**2 #formula for kinetic E
    
    r = np.sqrt(x**2 + y**2)  #position
    
    pot_energy = -G * m_s * m * (1 - beta) / r #potential energy
    
    return kinetic_energy , pot_energy


    
    