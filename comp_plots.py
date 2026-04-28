import matplotlib.pyplot as plt
import numpy as np
from orbital_elements import ecc_math
from config import t5 , t6 , t7

"""Plotting eccentricity from math and scaled expression"""
def ecc_comps(x , y , vx , vy , beta , t , type = None):
    """input: x (array), x position
              y (array), y position
              vx (array), x velocity
              vy (array), y velocity
              beta (array), beta values
              t (array), time values
              type (string), default: None, "numerical" if using numerical solved params
              else "perturbed", if using perturbed params"""
    
    theta = np.atan2(y , x) #angle
    theta_cont = np.unwrap(theta) / (2 * np.pi) #continuous angle

    # ecc_sc1 = eccentricity_sc(x , y , vx , vy , beta , t , type = type) #scaled eccentricity

    ecc_m , _ = ecc_math(x , y , t) #mathematical eccentricity

    ecc_m = ecc_m[::50]
    that = t[::50]

    orbit = np.floor(theta_cont / (2*np.pi)).astype(int) #iterating over orbits

    theta_cont = theta_cont[::50]
    orbit = orbit[::50]
    # ecc_sc1 = ecc_sc1[::50]
    
    # plt.plot(theta_cont , ecc_sc1 * 10**5 , label = "Scaled e")
    plt.plot(theta_cont , ecc_m , label = "Geometric e")
    
    plt.xlabel(r"$\hat{\theta}$ / $2\pi$")
    plt.ylabel(r"e value $\times 10^5$")
    plt.title("Scaled and geometric e as function of orbits")
    #plt.legend(bbox_to_anchor = (0.1 , 0.6))
    plt.legend()
    plt.show()


if __name__ == "__main__":
    rk = np.load("Files/rk45_t6_large_silicate_slowsw_realbeta.npz")
    x1 , y1 , vx1 , vy1 , m1 , b1 , t1 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b" , "t")]
    ecc_comps(x1 , y1 , vx1 , vy1 , b1 , t1)
    
    
    