import numpy as np

"""File calculating eccentricity based on scaled formula as well as using elliptical parameters"""

"""Scaled eccentricity formula"""
def ecc_scaled(num , particle_obj , pert = None):
    """input: num (tuple), consisting of numerically calculated x, y, vx, vy, b
              B (float), initial beta value
              pert (tuple), default:None, if perturbation evaluated pert contains
                                    perturbed r, omega and vr
        
       output: eccentricity (float)"""
    
    x , y , vx , vy , _ , b , _ = num

    theta = np.atan2(y , x) #numerical angle
    theta_cont = np.unwrap(theta) #continuous angle

    r = np.sqrt(x**2 + y**2) #numerical radial distance

    angvel_num = (-vx * np.sin(theta_cont) + vy * np.cos(theta_cont)) / r #numerical angular vel
    v_r = vx * np.cos(theta_cont) + vy * np.sin(theta_cont) #cartesian to radial numerical vel

    ecc_x = (-(1 - particle_obj.B) * r**2 * v_r * angvel_num / (1 - b * particle_obj.B) * (-y / r) - 
             (x / r) + r**3 * angvel_num**2 * (1 - particle_obj.B) / (1 - b * particle_obj.B) * (x / r))
    
    ecc_y = (-(1 - particle_obj.B) * r**2 * v_r * angvel_num / (1 - b * particle_obj.B) * (x / r) 
             - (y / r) + r**3 * angvel_num**2 * (1 - particle_obj.B) / (1 - b * particle_obj.B) * (y / r))    
    
    ecc_tot = np.sqrt(ecc_x**2 + ecc_y**2)

    #Iterating over orbits, one orbit is 2pi*n, n positive integer
    orbit_idx = np.floor(theta_cont / (2 * np.pi)).astype(int)
    change_idx = np.where(np.diff(orbit_idx) != 0)[0] + 1
    change_idx = np.insert(change_idx , 0 , 0)
    change_idx = np.append(change_idx, len(theta_cont) - 1)

    eccentricity = ecc_tot[change_idx] #one ecc sample each start of orbit
    orb = np.arange(1, len(change_idx) + 1) #updating orbits

    ecc_pert = None

    if pert is not None:
        r_per , omega_per , vr_per = pert #perturbed parameters

        ecc_theta = -(1 - particle_obj.B) / (1 - b * particle_obj.B) * r_per**2 * vr_per * omega_per #perturbed ecc theta dir
        ecc_r = (1 - particle_obj.B) / (1 - b * particle_obj.B) * r_per**3 * omega_per**2 - 1 #perturbed ecc radial dir

        eccpert_tot = np.sqrt(ecc_theta**2 + ecc_r**2) #total ecc
        ecc_pert = eccpert_tot[change_idx] #one ecc each start of orbit
    
    return eccentricity , ecc_pert , orb


"""Eccentricity from elliptical parameters"""
def ecc_calcs(x , y , r_pert = None):
    """input: x (array), numerical x position 
              y (array), numerical y position
              r_pert (array), perturbed radial position, default: None

       returns: eccenctricity, ecc_pert, orb (tuple), numerical and perturbed eccentricity and 
                number of orbits"""

    theta = np.atan2(y , x) #numerical angle
    theta_cont = np.unwrap(theta) #continuous angle

    r = np.sqrt(x**2 + y**2) #numerical radial distance
    
    start = 0
    k = 1
    eccentricity = [] #initiating numerical ecc list

    ecc_pert = [] #initiating perturbed ecc list
    orb = [] #initiating orb counts list

    #iterating through all angles
    for i in range(1 , len(theta_cont)):

        #for each full orbit, finding rmax and rmin, calculating numerical eccentricity
        if theta_cont[i] >= k * 2 * np.pi:
            segment = r[start : i]
            r_max = np.max(segment)
            r_min = np.min(segment)
            a = (r_max + r_min) / 2
            b = np.sqrt(r_max * r_min)
            eccentricity.append(np.sqrt(1 - b**2 / a**2))

            orb.append(k)

            #same calcs for perturbed equivalent  
            if r_pert is not None:
                segment_pert = r_pert[start : i]
                r_pertmax = np.max(segment_pert)
                r_pertmin = np.min(segment_pert)

                apert = (r_pertmax + r_pertmin) / 2
                bpert = np.sqrt(r_pertmax * r_pertmin)
                ecc_pert.append(np.sqrt(1 - bpert**2 / apert**2))

            start = i
            k += 1

    return eccentricity , ecc_pert , orb

    
    
    
    
    

    
    


