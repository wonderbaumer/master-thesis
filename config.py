import numpy as np
from scipy.constants import G
from constants import sputtering_yield , sw_flux , mA_S , rho , m_s , au , c , mA_C
from polar_to_cart import polar_to_cartesian

"""evaluating particle size against beta0"""
r_vals = np.linspace(0.00100 * 10**(-6) , 50 * 10**(-6) , 200)
r_betatest = np.linspace(0.00100 , 50 , 200)
m_range = 4 / 3 * np.pi * rho * r_vals**3 #masses corresponding to size range

"""calculates small parameter epsilon, using constants from mass calcs"""
def eps(m , material = "silicate" , sw = "slow" , species = "all"):
    """input: m (float), default M, mass of particle in kg
              sw (string), default: slow, options fast and CME
              species (string), default: all, else one of the elements H, He, C, O, N, Fe, Ne, Mg, Si, S

       returns: eps (float), epsilon parameter"""
    
    Ytot = sputtering_yield(material , sw , species) #total sputtering yield
    fsw = sw_flux(sw) #solar wind flux
    if material == "silicate":
        mA = mA_S
    
    elif material == "carbon": 
        mA = mA_C
    
    else: 
        raise ValueError("Material must be silicate or carbon")
    
    eps = fsw * Ytot * mA * np.pi * (3 / (4 * np.pi * rho))**(2 / 3) * m**(-1 / 3) * T #mass constant

    return eps

"""t hat combinations used, dt timestep, t_tot total simulation time"""
#1000 orbits
dt5 = 3.16e3 / 6834848 #original T
t_tot5 = 1000
t5 = (dt5 , t_tot5)

#10000 orbits
dt6 = 3.16e3 / 6834848
t_tot6 = 10000
t6 = (dt6 , t_tot6)

#20000 orbits
dt7 = 3.16e3 / 6834848
t_tot7 = 20000
t7 = (dt7 , t_tot7)


    