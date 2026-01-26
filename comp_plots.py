import matplotlib.pyplot as plt
import numpy as np
from orbital_elements import eccentricity_sc, ecc_math
from config import t5 , t6 , t7

def ecc_comps(x , y , vx , vy , beta , t):
    ecc_sc1 = eccentricity_sc(x , y , vx , vy , beta)
    orbits , math_arr , ecc_sc2 = ecc_math(ecc_sc1 , t , x , y)
    ecc_m = math_arr[: , 2]
    print(ecc_m)
    plt.plot(orbits , ecc_sc2 * 10**5 , label = "Scaled e")
    plt.plot(orbits , ecc_m * 10**5 , label = "Definition e")
    
    plt.xlabel("Full orbits")
    plt.ylabel(r"e value $\times 10^5$")

    plt.title("Scaled and definition e as function of full orbits")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    rk = np.load("Files/rk45_t6_masslossTrue_scaledeqs.npz")
    x1 , y1 , vx1 , vy1 , m1 , b1 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    ecc_comps(x1 , y1 , vx1 , vy1 , b1 , t6)
    #r = np.sqrt(x1**2+y1**2)
    #print(r[0:13582])