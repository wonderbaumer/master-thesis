import numpy as np
from scipy.interpolate import PchipInterpolator as pchip
import matplotlib.pyplot as plt
from config import B , eps , init_cart_scaled , delta , m_range , r_vals , r_init , M , V , c
from constants import sil_beta , car_beta , dat_to_arr , rho

sil_size , sil_betaval , _ = dat_to_arr(sil_beta)

size = sil_size * 10**(-6) / r_init
betaval = sil_betaval / B

interp = pchip(size , betaval)
b = interp(size)
# plt.plot(size , b)
# plt.xlabel("Particle size")
# plt.ylabel(r"$\hat{\beta}$ value")
# plt.title(r"$\hat{\beta}$ curve")
# plt.show()

def beta_real(par_size):

    b = interp(par_size)

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
def betahat(m):
    """input: m (float), scaled mass of particle

       returns: betahat(float), scaled betahat """

    r = m**(1/3)

    #b = m**(-1 / 3) #scaled betahat

    b = beta_real(r)

    return b

def betahat_car(m):
    """input: m (float), scaled mass of particle

       returns: betahat(float), scaled betahat for carbon """

    r = (3 * m / (4 * rho * np.pi))**(1 / 3)

    car_size , car_betaval , _ = dat_to_arr(car_beta)
    car_size = car_size * 10**(-6) #m
    
    #betahat = np.interp(r , car_size , car_betaval) #interpolating betahat for carbon
    log_size = np.log10(car_size)
    log_beta = np.log10(car_beta)

    #log_interp = interp1d(log_size, log_beta, kind='linear', fill_value='extrapolate')
    #betahat = 10**log_interp(np.log10(r))
    
    #return betahat

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

    r_hat_x = x / r
    r_hat_y = y / r
    v_dot_r = vx*r_hat_x + vy*r_hat_y

    ax = - betahat(m) * (vx + v_dot_r * r_hat_x) * B * V / (c * r**2 * (1 - B)) 
    ay = - betahat(m) * (vy + v_dot_r * r_hat_y) * B * V / (c * r**2 * (1 - B)) 

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
    
    