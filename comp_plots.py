import matplotlib.pyplot as plt
import numpy as np
from orbital_elements import eccentricity_sc, ecc_math
from config import t5 , t6 , t7 , c , V , B , au , T , t8

def ecc_comps(x , y , vx , vy , beta , t , type):
    theta = np.atan2(y , x) #angle
    theta_cont = np.unwrap(theta) / (2 * np.pi) #continuous angle

    ecc_sc1 = eccentricity_sc(x , y , vx , vy , beta , t , type = type)

    ecc_m , time_arr , ecc_sc1 = ecc_math(x , y , t , ecc_sc1)
    ecc_m = ecc_m[::50]
    

    that = time_arr[::50]
    orbit = np.floor(theta_cont / (2*np.pi)).astype(int)
    theta_cont = theta_cont[::50]
    orbit = orbit[::50]

    ecc_sc1 = ecc_sc1[::50]
    
    plt.plot(theta_cont , ecc_sc1 * 10**5 , label = "Scaled e")
    plt.plot(theta_cont , ecc_m * 10**5, label = "Geometric e")
    
    plt.xlabel(r"$\hat{\theta}$ / $2\pi$")
    plt.ylabel(r"e value $\times 10^5$")
    plt.title("Scaled and geometric e as function of orbits")
    #plt.legend(bbox_to_anchor = (0.1 , 0.6))
    plt.legend()
    plt.show()

def test_plot(x , y , vx , vy , m , b , t):
    dt , t_tot = t
    time = np.arange(0 , t_tot , dt)
    r = np.sqrt(x**2 + y**2)

    plt.plot(time / T , r / au , label = "r vs time")
    plt.show()




if __name__ == "__main__":
    rk = np.load("Files/rk45_t8_masslossFalse.npz")
    x1 , y1 , vx1 , vy1 , m1 , b1 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    #ecc_comps(x1 , y1 , vx1 , vy1 , b1 , t6 , "perturbed")
    test_plot(x1 , y1 , vx1 , vy1 , m1 , b1 , t8)
    
    
    