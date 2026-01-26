import numpy as np
import matplotlib.pyplot as plt
from config import B , t6 

"""func calculating scaled eccentricity using scaled/hatted values from numerical solver"""
def eccentricity_sc(x , y , vx , vy , beta):
    """input: x (array), x position 
              y (array), y position 
              vx (array), x velocity 
              vy (array), y velocity 
              beta (array), beta value 

       returns: ecc (arr), scaled eccentricity values"""
    
    r = np.sqrt(x**2 + y**2) #r hat

    theta_num = np.atan2(y , x) #theta hat
    theta_num = np.unwrap(theta_num) #continuous theta hat
    vr = vx * np.cos(theta_num) + vy * np.sin(theta_num) #cart to radial vel
    vtheta = -vx * np.sin(theta_num) + vy * np.cos(theta_num) #cart to theta vel
    omega = vtheta / r #angular velocity

    l_frac = 2 * r**4 * (1 - B)**2 * omega**2 / (1 - B * beta)**2 
    energies = 1 / 2 * vr**2 + 1 / 2 * r**2 * omega**2 - (1 - beta * B) / ((1 - B) * r)

    ecc_sq = 1 + l_frac * energies #eccentricity squared
    
    #setting extremely small negative numbers to zero
    small_val = 1e-14
    ecc_sq = np.where(ecc_sq<small_val , 0.0 , ecc_sq)

    ecc = np.sqrt(ecc_sq) #eccentricity

    return ecc

"""calculates eccentricity for each orbit based on mathematical definitions"""
def ecc_math(ecc_sc , t , x , y):
    """input: ecc_sc (array), scaled eccentricity values
              t (tuple), on the form dt, t_tot
              x (array), scaled x position 
              y (array), scaled y position
       returns: orbits (array), ecc_arr (array), ecc_scaled (array), 
       orbits, mathematical ecc and modified scaled ecc"""
    
    dt , t_tot = t #time unpacking
    that = np.arange(0 , t_tot , dt) #t hat

    orbit = round(2 * np.pi * len(that) / t_tot) #steps per 1 orbit
    r = np.sqrt(x**2 + y**2) #r hat
    orbits = np.arange(0 , len(that) - orbit + 1 , orbit) #orbits array, 1 element is 1 orbit

    N = int(len(r) / orbit)  #number of timesteps 

    ecc_arr = np.zeros((N , 3)) #math defs
    ecc_scaled = np.zeros((N , 1)) #scaled eqs
  
    row = 0

    #iterating over all orbits
    for i in orbits:
        if i+orbit > len(r): #running until orbits more than values for the orbit
            break

        r_orb = r[i:i+orbit] #r values for each orbit
        ecc_orb = ecc_sc[i:i+orbit] #scaled ecc values for each orbit

        r_max = np.max(r_orb) #largest radial distance
        r_min = np.min(r_orb) #smallest radial distance

        amajor = (r_max + r_min) / 2 #semi-major axis
        aminor = np.sqrt(r_max * r_min) #semi-minor axis
    
        ecce = (r_max - r_min) / (r_max + r_min)

        #semi-major, semi-minor and eccentricity into same list, different cols
        ecc_arr[row , 0] , ecc_arr[row , 1] , ecc_arr[row , 2] = r_min , r_max , ecce
        ecc_scaled[row] = ecc_orb[-1] #using final value for each orbit in scaled ecc as scaled ecc

        row +=1
    
    return orbits , ecc_arr , ecc_scaled

if __name__ == "__main__":
    rk = np.load("Files/rk45_t6_masslossTrue_scaledeqs.npz")
    x1 , y1 , vx1 , vy1 , m1 , b1 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    ecc = eccentricity_sc(x1 , y1 , vx1 , vy1 , b1)
    
    orbits , math , escaled = ecc_math(ecc , t6 , x1 , y1)
    print(math[: , 0] , math[: , 1])


