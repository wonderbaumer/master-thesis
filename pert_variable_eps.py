from dust_properties import dust_properties
from scipy.integrate import solve_ivp
from config import t6 , t7 , t8 , t9 , t10
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
    dbeta0_dt = beta0**2 / (3 * r0**2) #zeroth order perturbed dbeta0dt1

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

    return r0 , omega0 , beta0 , sol.t

if __name__== "__main__":
    par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "A")
    res = np.load("Files/rk45_t6_A_silicate_slowsw.npz")
    x , y , _ , _ , m , b , t  = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    t1 = par.epsilon * t
    func = rhs(par.epsilon * t , y0 , par.B , par.K)
    r0 , omega0 , beta0 , tpert = pert_motion(rhs , par.epsilon * t , y0 , par.B , par.K)
    
    plt.plot(t , r0)
    plt.plot(t , np.sqrt(x**2+y**2))
    plt.show()
    
    