import numpy as np
from scipy.constants import G
from constants import sputtering_yield , sw_flux , mA_S , rho , m_s , au , c , mA_C
from polar_to_cart import polar_to_cartesian

"""evaluating particle size against beta0"""
#r_vals = np.linspace(500e-9 , 10e-6 , 10) #size range over which to evaluate beta0
r_vals = np.linspace(1e-9 , 50e-6 , 200)
m_range = 4 / 3 * np.pi * rho * r_vals**3 #masses corresponding to size range
#M = 4 / 3 * np.pi * rho * (10.33226 * 10**(-6))**3 #mass corresponding to size of 10.33226 microns
#B = 0.0163 #beta value corresponding to size of 10.33226 microns
"""scaling parameters"""
B = 0.45931933916320633 #initial beta value
R = 1 * au #initial radial position
V = np.sqrt((G * m_s * (1 - B)) / R) #initial angular velocity, scaled formula
M = 1.30899694e-15 #initial particle mass in kg
T = round(np.sqrt(R**3 / (G * m_s * (1 - B)))) #initial period, scaled formula
K = 1

"""scaled initial parameters"""
rhat0 = R / R #initial scaled radial position
thetahat0 = 0 #initial scaled angular position
vrhat0 = 0  #initial scaled radial velocity
omegahat0 = V / V #initial scaled angular velocity
betahat0 = B / B #initial scaled beta
mhat0 = M / M #initial scaled mass

init_polar = np.array([R , 0 , 0 , V]) #initial cart coords
init_cart = polar_to_cartesian(init_polar) #initial cart coords

init_polar_scaled = np.array([rhat0 , thetahat0 , vrhat0 , omegahat0]) #initial scaled polar coords
init_cart_scaled = polar_to_cartesian(init_polar_scaled) #initial scaled cart coords
x , y , vx , vy = init_cart_scaled #unpacking init scaled cartesian coords

"""calculates small parameter epsilon, using constants from mass calcs"""
def eps(material = "silicate" , sw = "slow" , species = "all" , m = M):
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

delta = V / c

"""t hat combinations used, dt timestep, t_tot total simulation time"""
#1000 orbits
dt5 = 3.16e3 / T 
t_tot5 = 1000
t5 = (dt5 , t_tot5)

#10000 orbits
dt6 = 3.16e3 / T 
t_tot6 = 10000
t6 = (dt6 , t_tot6)

#20000 orbits
dt7 = 3.16e3 / T
t_tot7 = 20000
t7 = (dt7 , t_tot7)

if __name__ == "__main__":
    eps("carbon" , "slow")
    