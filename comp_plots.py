import matplotlib.pyplot as plt
import numpy as np
from orbital_elements import eccentricity_sc, ecc_math
from config import t5 , t6 , t7

def ecc_comps(x , y , vx , vy , beta , t):

    ecc_sc1 = eccentricity_sc(x , y , vx , vy , beta , t)

    ecc_m , time_arr , ecc_sc1 = ecc_math(x , y , t , ecc_sc1)
    ecc_m = ecc_m[::10]
    

    that = time_arr[::10]

    ecc_sc1 = ecc_sc1[::10]
    
    plt.plot(that , ecc_sc1 * 10**5 , label = "Scaled e")
    plt.plot(that , ecc_m * 10**5 , label = "Definition e")
    
    plt.xlabel("Number of orbits")
    plt.ylabel(r"e value $\times 10^5$")

    plt.title("Scaled and definition e as function of orbits")
    plt.legend(bbox_to_anchor = (0.1 , 0.6))
    #plt.legend()
    plt.show()

dt , t_tot = t7
t = np.arange(0 , t_tot , dt)
y = np.exp(t * 10**(-5))
plt.plot(t , 1-y) #huske å gange med Bbeta1/1-B for full effekt
plt.show()

if __name__ == "__main__":
    rk = np.load("Files/rk45_t6_masslossTrue_scaledeqs.npz")
    x1 , y1 , vx1 , vy1 , m1 , b1 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    