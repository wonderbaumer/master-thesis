import numpy as np
from config import B, M, V , R , T

"""function that calculates the scaled total energy of the particle in orbit, cartesian coordinates"""
def tot_energy(x , y , vx , vy ,  m , beta):
    """input: x (float), scaled x position
              y (float), scaled y position
              vx (float), scaled x velocity
              vy (float), scaled y velocity
              m (float), scaled mass of particle
              beta (float), betahat value
              
        returns: kinetic_energy, pot_energy (tuple), the scaled energies of the particle
        """
        
    v = np.sqrt(vx**2 + vy**2) #speed
    theta_num = np.atan2(y , x) #thetahat
    theta_num = np.unwrap(theta_num) #avoiding discontinuities
    r = np.sqrt(x**2 + y**2)  #position
    vr = (x*vx + y*vy)/r
    vtheta = (x*vy - y*vx) / r
    
    kinetic_energy = 1 / 2 * m * (vr**2 + vtheta**2)
    #kinetic_energy = 1 / 2 * m * v**2 #formula corresponding to scaled kinetic E
    
    r = np.sqrt(x**2 + y**2)  #position
    
    pot_energy = -(1 - B * beta) * m / ((1 - B) * r) #scaled potential energy
    
    return kinetic_energy , pot_energy
