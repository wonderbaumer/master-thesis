from dust_properties import dust_properties
from scipy.integrate import solve_ivp
from config import t6 , t7 , t8 , t9 , t10 , yr , lower_m_sil , lower_m_car
import numpy as np
import matplotlib.pyplot as plt


"""This code solves the zeroth order perturbed expressions (multiscale) numerically"""

y0 = [1.0 , 1.0] #Initial conditions beta0 and C0

"""Differential equations to solve"""
def rhs(t , y0 , B , K):
    """Inputs: t (array), simulation time
               y0 (list), initial conditions
               B (float), initial beta value
               K (float), drag-to-mass loss ratio

       Returns: diff_eqs (list), expressions to solve for beta0 and C0"""
    
    beta0 = y0[0]
    C0 = y0[1]

    m0 = 1.0
    beta0 = m0**(-1 / 3)

    r0 = C0**2 * (1 - B) / (1 - B * beta0) #zeroth order perturbed sol r0
    dbeta0_dt = beta0**2 / 3 * r0**(-2) #zeroth order perturbed dbeta0dt1

    dC0_dt = -B * K * beta0 * (1 - beta0 * B)**2 / (1 - B)**3 * C0**(-3) #zeroth order perturbed dC0dt1
    diff_eqs = [dbeta0_dt , dC0_dt]

    return diff_eqs

"""Stopping condition if size outside interpolation range"""
def r_out_of_range_event(t , y , B , K):
    beta0 , C0 = y

    m0 = 1.0
    lower_lim = 1e-23

    mass_denom = 1 - B * beta0
    
    return mass_denom - 1e-2

r_out_of_range_event.terminal = True
r_out_of_range_event.direction = -1

"""Stopping condition if distance <0.1"""
def orbital_radius_event(t , y , B , K):
    beta0 , C0 = y

    r0 = C0**2 * (1 - B) / (1 - B * beta0)

    return r0 - 1e-2

orbital_radius_event.terminal = True
orbital_radius_event.direction = -1


"""Numerically solving the differential equations"""
def pert_motion(fun , t , y0 , B , K):
    """Inputs: fun (function), containing the differential equations to solve
               t (array), simulation time
               y0 (list), initial values
               B (float), initial beta value
               K (float), drag-to-mass loss ratio
            
        Returns: r0, omega0, beta0, t (tuple), solved parameters"""
    
    sol = solve_ivp(fun , (t[0] , t[-1]) , y0 , t_eval = t , args = (B , K) ,  method = "RK45" 
                    , rtol = 1e-6 , atol = 1e-9 , events = [r_out_of_range_event , orbital_radius_event])

    return sol

def arr_variables(sol , B):
    beta0 , C0 = sol.y

    r0 = C0**2 * (1 - B) / (1 - B * beta0)
 
    omega0 = C0**(-3) * ((1 - B) / (1 - B * beta0))**(-2)
    mpert = beta0**(-3)

    dbeta0_dt1 = beta0**2 / 3 * r0**(-2)

    K_exp = dbeta0_dt1 * (1 - B) / (2 * beta0 * (1 - beta0 * B))

    return r0 , omega0 , beta0 , sol.t , mpert , K_exp

if __name__== "__main__":
    material = "carbon"
    par = dust_properties(material , "CME" , init_dist = 1.0 , size = "C")
    # res = np.load("Files/rk45_t8_A_silicate_fastsw.npz")
    # x , y , _ , _ , m , b , t  = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    dt , t_tot = t9
    t = np.arange(0 , t_tot , dt)
    
    t1 = par.epsilon * t
    # func = rhs(par.epsilon * t , y0 , par.B , par.K)
    pert = pert_motion(rhs , par.epsilon * t , y0 , par.B , par.K)
    
    stopping_reason = None
    stopping_events = ["particle outside interpolation range" , "particle impacted Sun"]
    for i , te in enumerate(pert.t_events):
                if len(te) > 0:
                    stopping_reason = stopping_events[i]
                    break

    print(stopping_reason)

    r0 , omega0 , beta0 , tpert , mpert , Kpert = arr_variables(pert , par.B)
    
    m_phys = mpert * par.m0
    # print(tpert[-1] / par.epsilon * par.T / yr , r0[-1] , m_phys[-1] , beta0[-1])
    
    # plt.plot(tpert / par.epsilon , beta0)
    # plt.ylim(0 , 2.0)
    # plt.show()
    
    print(Kpert)
    
    