import numpy as np
from config import sil_betaval , car_betaval_bound , init_vals , material_files , m_s , R , size_to_mass , car_size_bound , sil_size
from scipy.constants import c , G
from forces_scaled import inter_func
from dust_properties import dust_properties
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator as pchip
from numpy.lib.stride_tricks import sliding_window_view
from config import t6 , t5
import scipy.integrate as it
import matplotlib.pyplot as plt

def perturbed_functions(t , B , K):

    def rhs(t , y):
        beta0 = y[0]
        C0 = y[1]

        m0 = 1.0
        beta0 = m0**(-1 / 3)

        r0 = C0**2 * (1 - B) / (1 - B * beta0)
        dbeta0_dt = beta0**2 / (3 * r0**2)

        dC0_dt = -B * K * beta0 * (1 - beta0 * B)**2 / (1 - B)**3 * C0**(-3)

        return [dbeta0_dt , dC0_dt]

    y0 = [1.0 , 1.0]

    sol = solve_ivp(rhs , (t[0] , t[-1]), y0, t_eval=t, method="RK45", rtol=1e-9 , atol=1e-12 )

    beta0 = sol.y[0]
    C0 = sol.y[1]


    r0 = C0**2 * (1 - B) / (1 - B * beta0)
 
    omega0 = C0**(-3) * ((1 - B) / (1 - B * beta0))**(-2)
    

    return r0 , omega0 , beta0 , sol.t

if __name__== "__main__":
    par = dust_properties("silicate" , "slow" , "large")
    res = np.load("Files/rk45_t6_large_silicate_slowsw.npz")
    x , y , _ , _ , m , b , t  = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]

    
    rad , omega0 , beta0 , time = perturbed_functions(par.epsilon0 * t , par.B , par.K)
   
    
    
    