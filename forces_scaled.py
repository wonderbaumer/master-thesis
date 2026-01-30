import numpy as np
from config import B , eps , init_cart_scaled

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

    betahat = m**(-1 / 3) #scaled betahat
    
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

"""function that calculates total acceleration based on pressure radiation force and gravity"""
def tot_acc(x , y , m):
    """input: x (float), scaled x position
              y (float), scaled y position
              m (float), scaled mass of particle

        returns: ax , ay (array), scaled acceleration in x and y dir"""
    
    px , py = pressure_radial(x , y , m) #decomposing pressure radiation force
    gx , gy = gravity(x , y) #decomposing gravitational force

    ax = px + gx  #total acceleration in x dir
    ay = py + gy  #total acceleration in y dir 
    
    return ax , ay

if __name__ == "__main__":
    x , y , vx , vy = init_cart_scaled

    a = sputtering(1e-15 , sw = "fast" , species = "H")
    print(a)
    