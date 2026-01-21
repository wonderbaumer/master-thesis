import numpy as np
import matplotlib.pyplot as plt
from config import B , t6
from energy import tot_energy

def eccentricity(x , y , vx , vy , m , beta):
    k , u = tot_energy(x , y , vx , vy ,  m , beta)
    e = k + u

    r = np.sqrt(x**2 + y**2)
    v = np.sqrt(vx**2 + vy**2)
    omega = r / v
    l = omega * m * r**2

    alpha =-(1 - B * beta) / (1 - B)

    #ecc = np.sqrt(1 + 2 * e * l**2 / (m * alpha**2))
    ecc = 1 + 2 * e * l**2 / (m * alpha**2)

    return ecc

if __name__ == "__main__":
    rk = np.load("Files/rk45_t6_masslossTrue_scaledeqs.npz")
    x1 , y1 , vx1 , vy1 , m1 , b1 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    ecc = eccentricity(x1 , y1 , vx1 , vy1 , m1 , b1)

    dt , t_tot = t6
    that = np.arange(0 , t_tot , dt)

    orbit = round(2 * np.pi * len(that) / t_tot)

    theta = np.atan2(y1 , x1) #rk45 thetahat
    theta = np.unwrap(theta)
    

    orbits = np.arange(orbit , len(that) , orbit)
    """
    for i in orbits:

        theta = theta[i:i+1]

        x = x1[i:i+1]
        y = y1[i:i+1]
        theta = theta[i:i+1]

        i_arr = 
        #a = x1[10:] / np.cos(theta[10:])
        #b = y1[10:] / np.sin(theta[10:])

        #a_arr = np.array(a)
        #b_arr = np.array(b)
    """


    

