import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from config import B , eps , init_cart_scaled , delta , m_range , r_vals
from constants import sil_beta , car_beta , dat_to_arr , rho

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
def betahat(m):
    """input: m (float), scaled mass of particle

       returns: betahat(float), scaled betahat """

    #betahat = m**(-1 / 3) #scaled betahat

    r = (3 * m / (4 * rho * np.pi))**(1 / 3)

    sil_size , sil_betaval , _ = dat_to_arr(sil_beta)
    sil_size = sil_size * 10**(-6) #m
    betahat = np.interp(r , sil_size , sil_betaval) #interpolating betahat for silicate

    return betahat

def betahat_car(m):
    """input: m (float), scaled mass of particle

       returns: betahat(float), scaled betahat for carbon """

    r = (3 * m / (4 * rho * np.pi))**(1 / 3)

    car_size , car_betaval , _ = dat_to_arr(car_beta)
    car_size = car_size * 10**(-6) #m
    
    #betahat = np.interp(r , car_size , car_betaval) #interpolating betahat for carbon
    log_size = np.log10(car_size)
    log_beta = np.log10(car_beta)

    log_interp = interp1d(log_size, log_beta, kind='linear', fill_value='extrapolate')
    betahat = 10**log_interp(np.log10(r))
    
    return betahat

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

    b_term = betahat(m) * B / ((1 - B) * r**4) #beta term in expression 
    xvel_terms = -2 * vx * (x**2 + y**2) - x * y * vy #x velocity terms
    yvel_terms = -x * y * vx - 2 * vy * (y**2 + x**2) #y velocity terms

    ax = b_term * delta * xvel_terms #acceleration in x direction
    ay = b_term * delta * yvel_terms #acceleration in y direction

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
    
    plt.xscale("log")
    plt.yscale("log")
    plt.plot(r_vals , betahat(m_range))
    plt.show()

   
    