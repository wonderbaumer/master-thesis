import numpy as np
from config import B

"""function that calculates the scaled total energy of the particle in orbit, cartesian coordinates"""
def tot_energy(x , y , vx , vy , m , beta):
    """input: x (float), scaled x position
              y (float), scaled y position
              vx (float), scaled x velocity
              vy (float), scaled y velocity
              m (float), scaled mass of particle
              beta (float), betahat value
              
        returns: kinetic_energy, pot_energy (tuple), the scaled energies of the particle
        """
        
    v = np.sqrt(vx**2 + vy**2) #speed
    
    kinetic_energy = 1 / 2 * m * v**2 #formula corresponding to scaled kinetic E
    
    r = np.sqrt(x**2 + y**2)  #position
    
    pot_energy = -(1 - B * beta) / ((1 - B) * r) #scaled potential energy
    
    return kinetic_energy , pot_energy
