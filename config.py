from constants import *
from polar_to_cart import *

"""initial values"""
m_par = 1.30899694e-15 #mass particle in kg
r_par_init = 500e-9 #radius particle in m
beta0 = 0.45931933916320633 #initial mass and radius from initial pressure radiation/gravity
vtheta0 = 2.19013101e+04 #initial angular vel in ms^-1, only used for testing purposes

r0 = 1.0 * au #initial radial dist in units of AU
theta0 = 0 #initial angular position in rad
v0r = 0 #initial radial vel in ms^-1
init_polar = np.array([r0 , theta0 , v0r]) #initial values array

init_polar_test = np.array([r0 , theta0 , v0r , vtheta0]) #initial values array, only used for testing
init_cartesian = polar_to_cartesian(init_polar_test) #initial values to cartesian, only used for testing

"""evaluating particle size against beta0"""
r_vals = np.linspace(500e-9 , 10e-6 , 10) #size range over which to evaluate beta0
m_range = 4 / 3 * np.pi * rho * r_vals**3 #masses corresponding to size range

"""scaling parameters"""
R = au #initial radial position
V = np.sqrt((G * m_s * (1 - beta0)) / R) #initial angular velocity, scaled formula
T = np.sqrt(R**3 / (G * m_s * (1 - beta0))) #initial period, scaled formula

rhat0 = r0 / R #initial scaled radial position
vtheta0scaled = vtheta0 / V #initial scaled angular velocity

"""calculates small parameter epsilon, using constants from mass calcs"""
def eps(m = m_par):
    """input: m (float), default m_par, mass of particle in kg

       returns: eps (float), epsilon parameter"""
    
    eps = fsw * Ytot * mA * np.pi * (3 / (4 * np.pi * rho))**(2 / 3) * m**(-1 / 3) * T #mass constant

    return eps

orb_period = T * 2 * np.pi #orbital period in s
orb_per_yr = orb_period / yr #orbital period in years

"""time combinations used"""
#1 orbit
dt1 = 3.16e1
t_tot1 = orb_period
t1 = (dt1 , t_tot1)

#10 orbits
dt2 = 3.16e2
t_tot2 = 10 * orb_period
t2 = (dt2 , t_tot2)

#100 orbits
dt3 = 3.16e3
t_tot3 = 100 * orb_period
t3 = (dt3 , t_tot3)

#500 orbits
dt4 = 3.16e3
t_tot4 = 500 * orb_period
t4 = (dt4 , t_tot4)

if __name__ == "__main__":
    print(orb_period / yr)