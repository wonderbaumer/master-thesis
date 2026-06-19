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

"""Numerically solving the differential equations"""
def pert_motion(fun , t , y0 , B , K):
    """Inputs: fun (function), containing the differential equations to solve
               t (array), simulation time
               y0 (list), initial values
               B (float), initial beta value
               K (float), drag-to-mass loss ratio
            
        Returns: r0, omega0, beta0, t (tuple), solved parameters"""
    
    sol = solve_ivp(fun , (t[0] , t[-1]) , y0 , t_eval = t , args = (B , K) ,  method = "RK45" 
                    , rtol = 1e-9 , atol = 1e-12)

    beta0 = sol.y[0]
    C0 = sol.y[1]

    r0 = C0**2 * (1 - B) / (1 - B * beta0)
 
    omega0 = C0**(-3) * ((1 - B) / (1 - B * beta0))**(-2)
    mpert = beta0**(-3)

    return r0 , omega0 , beta0 , sol.t , mpert

if __name__== "__main__":
    material = "silicate"
    par = dust_properties(material , "fast" , init_dist = 1.0 , size = "D")
    res = np.load("Files/rk45_t8_A_silicate_fastsw.npz")
    # x , y , _ , _ , m , b , t  = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    dt , t_tot = t10
    t = np.arange(0 , t_tot , dt)
    
    t1 = par.epsilon * t
    func = rhs(par.epsilon * t , y0 , par.B , par.K)
    r0 , omega0 , beta0 , tpert , mpert = pert_motion(rhs , par.epsilon * t , y0 , par.B , par.K)

    # lower_m_mat = lower_m_sil if material == "silicate" else lower_m_car
    # m_phys = mpert * par.m0
    # idx = []
    # for i in range(0 , len(r0)-1):
    #     if r0[i] > 0.1:
    #         idx.append(i)

    # idx = idx[-1]
    # tp1 = tpert[idx] 
    # idxm = []
    # for i in range(0 , len(m_phys) - 1):
    #     if m_phys[i] > lower_m_mat:
    #         idxm.append(i)

    # idxm = idxm[-1] 

    # tp2 = tpert[idxm]
    # r = r0[idx]
    
    # m_phys = m_phys[m_phys>lower_m_mat]
    # print(tp1 / par.epsilon * par.T / yr , tp2 / par.epsilon * par.T / yr , r)
    
    #  , r , m_phys
    
    
    