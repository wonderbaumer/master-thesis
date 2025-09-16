import numpy as np
from forces import gravity
from polar_to_cart import polar_to_cartesian
from leapfrog import leapfrog_algorithm
from scipy.constants import *

"""scipy constants used:
    G: gravitational constant: 6.6743e-11 m^3kg^-1s^-2"""

m_s = 1.98847e30  #mass of sun, in kg

"""function that calculates the total energy of the particle in orbit,
uses cartesian coordinates, kinetic energy per unit mass"""
def tot_energy(x , y , vx , vy):
    """input: x (float), x position
              y (float), y position
              vx (float), x velocity
              vy (float), y velocity
              
        returns: kinetic_energy, pot_energy (float), energies
        """
        
    v = np.sqrt(vx**2 + vy**2) #speed
    
    kinetic_energy = 1 / 2 * v**2 #formula for kinetic E
    
    r = np.sqrt(x**2 + y**2) #position
    
    pot_energy = -G * m_s / r #potential energy
    
    return kinetic_energy , pot_energy