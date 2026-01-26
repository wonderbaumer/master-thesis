import numpy as np
import matplotlib.pyplot as plt
from config import B , t6 , T
from energy import tot_energy

def eccentricity(x , y , vx , vy , m , beta):
    r = np.sqrt(x**2 + y**2)

    theta_num = np.atan2(y , x) #thetahat
    theta_num = np.unwrap(theta_num)
    vr = vx * np.cos(theta_num) + vy * np.sin(theta_num)
    vtheta = -vx * np.sin(theta_num) + vy * np.cos(theta_num)
    omega = vtheta / r

    l_frac = 2 * r**4 * (1 - B)**2 * omega**2 / (1 - B * beta)**2
    energies = 1 / 2 * vr**2 + 1 / 2 * r**2 * omega**2 - (1 - beta * B) / ((1 - B) * r)

    ecc_sq = 1 + l_frac * energies

    small_val = 1e-14
    ecc_sq = np.where(ecc_sq<small_val , 0.0 , ecc_sq)

    ecc = np.sqrt(ecc_sq)

    return ecc

if __name__ == "__main__":
    rk = np.load("Files/rk45_t6_masslossTrue_scaledeqs.npz")
    x1 , y1 , vx1 , vy1 , m1 , b1 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    ecc = eccentricity(x1 , y1 , vx1 , vy1 , m1 , b1)
    
    dt , t_tot = t6
    that = np.arange(0 , t_tot , dt)

    #plt.plot(that , ecc)
    #plt.show()


    orbit = round(2 * np.pi * len(that) / t_tot) #steps per 1 orbit
    
    r = np.sqrt(x1**2 + y1**2)
    

    orbits = np.arange(0 , len(that) - orbit + 1 , orbit)

    N = int(len(r) / orbit)  #number of timesteps
    ecc_arr = np.zeros((N , 3)) #math defs
    ecc_scaled = np.zeros((N , 1)) #scaled eqs
  
    row = 0

    for i in orbits:
        if i+orbit > len(r):
            break

        r_orb = r[i:i+orbit]
        ecc_orb = ecc[i:i+orbit]

        r_max = np.max(r_orb)
        r_min = np.min(r_orb)

        amajor = (r_max + r_min) / 2
        aminor = np.sqrt(r_max * r_min)
    
        ecce = np.sqrt(1 - aminor**2 / amajor**2)

        ecc_arr[row , 0] , ecc_arr[row , 1] , ecc_arr[row , 2] = amajor , aminor , ecce
        ecc_scaled[row] = ecc_orb[-1]

        row +=1

#plt.plot(orbits , ecc_arr[: , 2] , label = "math defs")
#plt.plot(orbits, ecc_scaled , label = "scaled eqs")
#plt.legend()
#plt.show()
print(ecc_arr[: , 2] , ecc_scaled)   



    

