import numpy as np
import matplotlib.pyplot as plt
from config import B , t6 
from pert_functions import vrhat_pert , thetahat_pert , omegahat_pert , rhat_pert , betahat_pert

"""func calculating scaled eccentricity using scaled/hatted values from numerical solver"""
def eccentricity_sc(x , y , vx , vy , beta , t):
    """input: x (array), x position 
              y (array), y position 
              vx (array), x velocity 
              vy (array), y velocity 
              beta (array), beta value 

       returns: ecc (arr), scaled eccentricity values"""
    
    dt , t_tot = t
    that = np.arange(0 , t_tot , dt)

    #Numerical
    r = np.sqrt(x**2 + y**2) #r hat
    theta_num = np.atan2(y , x) #theta hat
    theta_num = np.unwrap(theta_num) #continuous theta hat
    vr = vx * np.cos(theta_num) + vy * np.sin(theta_num) #cart to radial vel
    vtheta = -vx * np.sin(theta_num) + vy * np.cos(theta_num) #cart to theta vel
    omega = vtheta / r #angular velocity

    #Perturbed
    r_p = rhat_pert(that)
    omega_p = omegahat_pert(that)
    vr_p = vrhat_pert(that)
    beta_p = betahat_pert(that)


    #l_frac = 2 * r**4 * (1 - B)**2 * omega**2 / (1 - B * beta)**2 #not pert
    #energies = 1 / 2 * vr**2 + 1 / 2 * r**2 * omega**2 - (1 - beta * B) / ((1 - B) * r) #not pert

    l_frac = 2 * r_p**4 * (1 - B)**2 * omega_p**2 / (1 - B * beta_p)**2
    energies = 1 / 2 * vr_p**2 + 1 / 2 * r_p**2 * omega_p**2 - (1 - beta_p * B) / ((1 - B) * r_p)

    ecc_sq = 1 + l_frac * energies #eccentricity squared
    
    #setting extremely small negative numbers to zero
    small_val = 1e-14
    ecc_sq = np.where(ecc_sq<small_val , 0.0 , ecc_sq)

    ecc = np.sqrt(ecc_sq) #eccentricity

    return ecc

"""calculates eccentricity for each orbit based on mathematical definitions"""
def ecc_math(x , y , t , ecc_sc):
    """input: t (tuple), on the form dt, t_tot
              x (array), scaled x position 
              y (array), scaled y position

       returns: eccearr_filled (array), eccentricity array corresponding to time values"""

    theta = np.atan2(y , x) #angle
    theta_cont = np.unwrap(theta) #continuous angle

    r = np.sqrt(x**2 + y**2) #radial distance

    dt , t_tot = t
    that = np.arange(0 , t_tot , dt)
  
    orbit_index = np.floor(theta_cont / (2 * np.pi)) #defining orbit in terms of angle
    full_orbit = np.where(np.diff(orbit_index) > 0)[0] + 1 #indexing orbits 

    orb_start = np.concatenate(([0] , full_orbit)) #start point of each orbit
    orb_end = np.concatenate((full_orbit , [len(r)])) #end point of each orbit
    counts = orb_end - orb_start #number of steps per orbit
    
    r_max = np.array([r[start : end].max() for start , end in zip(orb_start , orb_end)]) #max r val per orbit
    r_min = np.array([r[start : end].min() for start , end in zip(orb_start , orb_end)]) #min r val per orbit
    amajor = (r_max + r_min) / 2 #semi-major axis
    aminor = np.sqrt(r_max * r_min) #semi-minor axis

    ecce = np.sqrt(1 - (aminor**2 / amajor**2)) #eccentricity 

    ecce_arr = np.stack([aminor , amajor , ecce] , axis = 1) #array of aminor, amajor, eccentricity
    eccearr_filled = np.repeat(ecce_arr , counts , axis = 0) #filling remainder steps with equal values of eccentricity for comparing with that and other ecc

    ecce_filled = eccearr_filled[: , 2]

    return ecce_filled , that , ecc_sc


if __name__ == "__main__":
    rk = np.load("Files/rk45_t6_masslossTrue_scaledeqs.npz")
    x1 , y1 , vx1 , vy1 , m1 , b1 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    ecc = eccentricity_sc(x1 , y1 , vx1 , vy1 , b1 , t6)
    
    
    
    

    
    


