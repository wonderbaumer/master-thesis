import numpy as np
import matplotlib.pyplot as plt
from config import B , eps , init_cart_scaled , delta , m_range , r_vals , r1 , r2 , r3 , M , V , c
from scipy.interpolate import PchipInterpolator as pchip

def inter_func(self , bval_file):
        size , betaval , _ = dat_to_arr(bval_file)
        size *= 1e-6 / self.size

        return pchip(size , betaval)


def beta_real(par_size , material , b_val):
        if material == "Silicate":
            b = silicate_interp(par_size) / b_val
    
        else:
            b = carbon_interp(par_size) / b_val

        return b

"""calculates acceleration of the particle in x and y direction
    based on scaled gravitational force between particle and Sun"""
def gravity(x , y):
    """input: x (float), scaled x position
              y (float), scaled y position
              
       returns: acc_x, acc_y (float), scaled acceleration in x and y direction"""
       
    r = np.sqrt(x**2 + y**2) #radial distance
    acc_x = - x / ((1 - B) * r**3) #scaled acceleration in x direction
    acc_y = - y / ((1 - B) * r**3) #scaledacceleration in y direction

    return acc_x , acc_y 

"""calculates mass change from sputtering based on scaled parameters"""
def sputtering(m , epsilon):
    """input: m (float), scaled mass
        
       return: dmdt (float), mass change as function of time"""

    #dmdt = - eps(sw , species) * m**(2 / 3) #mass change with time
    dmdt = - epsilon * m**(2 / 3)
    
    return dmdt

"""function that calculates betahat, based on scaled equations"""
def betahat(m , particle_obj):
    """input: m (float), scaled mass of particle

       returns: betahat(float), scaled betahat """

    r = m**(1/3)

    #b = m**(-1 / 3) #scaled betahat

    #b = beta_real(r , material)
    b = particle_obj.beta_real(r)

    return b


"""function that calculates the radial component of the pressure radiation force, 
based on scaled equations"""
def pressure_radial(x , y , m):
    """input: x (float), scaled x position
              y (float), scaled y position
              m (float), scaled mass of particle

        returns: ax , ay (tuple), acceleration in x and y dir for pressure radiation force"""
    
    r = np.sqrt(x**2 + y**2) #scaled radial distance

    ax = x * betahat(m) * B / ((1 - B) * r**3) #scaled acceleration in x dir
    ay = y * betahat(m) * B / ((1 - B) * r**3) #scaled acceleration in y dir

    return ax , ay

"""calculates Poynting-Robertson drag based on scaled parameters"""
def pr_drag(x , y , vx , vy , m):
    """input: x (float), scaled x position
              y (float), scaled y position
              vx (float), scaled x velocity
              vy (float), scaled y velocity
              m (float), scaled particle mass

        returns ax , ay (tuple), acceleration in x and y direction"""
    
    r = np.sqrt(x**2 + y**2) #scaled radial distance

    theta = np.atan2(y , x)

    A = -betahat(m) * B * V / ((1 - B) * r**3 * c)
    x_dir = 2 * np.cos(theta) * (x * vx + y * vy) - np.sin(theta) * (x * vy - y * vx)
    y_dir = 2 * np.sin(theta) * (x * vx + y * vy) + np.cos(theta) * (x * vy - y * vx)

    ax = A * x_dir
    ay = A * y_dir

    return ax , ay

"""function that calculates total acceleration based on pressure radiation force, gravity and Poynting-Robertson drag"""
def tot_acc(x , y , vx , vy , m):
    """input: x (float), scaled x position
              y (float), scaled y position
              vx (float), scaled x velocity
              vy (float), scaled y velocity
              m (float), scaled mass of particle

        returns: ax , ay (array), scaled acceleration in x and y dir"""
    
    px , py = pressure_radial(x , y , m) #decomposing pressure radiation force
    gx , gy = gravity(x , y) #decomposing gravitational force
    prx , pry = pr_drag(x , y , vx , vy , m) #decomposing Poynting-Robertson force

    ax = px + gx + prx  #total acceleration in x dir
    ay = py + gy + pry  #total acceleration in y dir 
    
    return ax , ay

if __name__ == "__main__":
    x , y , vx , vy = init_cart_scaled
    
    