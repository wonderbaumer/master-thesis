import numpy as np
import matplotlib.pyplot as plt
from pert_functions import perturbed_functions
from dust_properties import dust_properties
from config import R

def ecc_scaled(num , B , pert = None):
    "input: pert (tuple), rpert,omegapert,vrpert"
    
    x , y , vx , vy , _ , b , _ = num

    theta = np.atan2(y , x) #angle
    theta_cont = np.unwrap(theta) #continuous angle

    r = np.sqrt(x**2 + y**2)
    start = 0
    k = 1

    # eccentricity = []
    # ecc_pert = []
    # orb = []

    angvel_num = (-vx * np.sin(theta_cont) + vy * np.cos(theta_cont)) / r
    v_r = vx * np.cos(theta_cont) + vy * np.sin(theta_cont) #cartesian to radial vel

    ecc_x = -(1 - B) * r**2 * v_r * angvel_num / (1 - b * B) * (-y / r) - (x / r) + r**3 * angvel_num**2 * (1 - B) / (1 - b * B) * (x / r)
    ecc_y = -(1 - B) * r**2 * v_r * angvel_num / (1 - b * B) * (x / r) - (y / r) + r**3 * angvel_num**2 * (1 - B) / (1 - b * B) * (y / r)
    ecc_tot = np.sqrt(ecc_x**2 + ecc_y**2)

    orbit_idx = np.floor(theta_cont / (2 * np.pi)).astype(int)
    change_idx = np.where(np.diff(orbit_idx) != 0)[0] + 1
    change_idx = np.insert(change_idx , 0 , 0)
    change_idx = np.append(change_idx, len(theta_cont) - 1)

    eccentricity = ecc_tot[change_idx]
    orb = np.arange(1, len(change_idx) + 1)

    ecc_pert = None

    if pert is not None:
        r_per, omega_per, vr_per = pert

        ecc_theta = -(1 - B) / (1 - b * B) * r_per**2 * vr_per * omega_per
        ecc_r = (1 - B) / (1 - b * B) * r_per**3 * omega_per**2 - 1

        eccpert_tot = np.sqrt(ecc_theta**2 + ecc_r**2)
        ecc_pert = eccpert_tot[change_idx]
    
    return eccentricity , ecc_pert , orb


"""calculates eccentricity for each orbit based on mathematical definitions"""
def ecc_calcs(x , y , r_pert = None):
    """input: t (tuple), on the form dt, t_tot
              x (array), scaled x position 
              y (array), scaled y position

       returns: eccearr_filled (array), eccentricity array corresponding to time values"""

    theta = np.atan2(y , x) #angle
    theta_cont = np.unwrap(theta) #continuous angle

    r = np.sqrt(x**2 + y**2) #radial distance
    
    start = 0
    k = 1
    eccentricity = []

    ecc_pert = []
    orb = []


    for i in range(1 , len(theta_cont)):

        if theta_cont[i] >= k * 2 * np.pi:

            segment = r[start : i]
            r_max = np.max(segment)
            r_min = np.min(segment)
            # eccentricity.append((r_max - r_min) / (r_max + r_min))
            a = (r_max + r_min) / 2
            b = np.sqrt(r_max * r_min)
            eccentricity.append(np.sqrt(1 - b**2 / a**2))

            orb.append(k)
                
            if r_pert is not None:
                segment_pert = r_pert[start : i]
                r_pertmax = np.max(segment_pert)
                r_pertmin = np.min(segment_pert)
                # ecc_pert.append((r_pertmax - r_pertmin) / (r_pertmax + r_pertmin))

                apert = (r_pertmax + r_pertmin) / 2
                bpert = np.sqrt(r_pertmax * r_pertmin)
                ecc_pert.append(np.sqrt(1 - bpert**2 / apert**2))

            start = i
            k += 1

    return eccentricity , ecc_pert , orb

if __name__ == "__main__":
    par = dust_properties("silicate" , "slow" , "large")
    file_path = "Files/rk45_t6_large_silicate_slowsw.npz"
    rk = np.load(file_path)
    x1 , y1 , vx1 , vy1 , m1 , b , t = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b" , "t")]
    p = perturbed_functions(par , t , b , find_k = False)
    c0 = p.C0(p.K)
    om , om0 , _ = p.omega(p.K)
    r_pert , r0 , r1 = p.rad(p.K)
    thetaval = p.theta(p.K)
    vrpert = p.vr(p.K)

    ecc , orb , ecc_pert = ecc_scaled((x1 , y1 , vx1 , vy1 , m1 , b , t) , par.B , (r_pert , om , vrpert))
    
    
    plt.plot(orb , ecc)
    plt.plot(orb , ecc_pert)
    plt.show()
    
    
    
    

    
    


