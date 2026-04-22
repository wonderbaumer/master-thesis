import numpy as np
import pandas as pd
from polar_to_cart import polar_to_cartesian

"""radiation pressure"""
rho_s = 3500 #kgm^-3 silicate
rho_c = 1800 #kgm^-3 carbon
S_s = 1361 #solar constant, in Wm^-2
au = 149597871e3 #one astronomical unit, AU, in m
q_pr = 1 #radiation pressure coefficient, unitless
m_s = 1.98847e30  #mass of sun, in kg
#c = 299792458 #speed of light, in ms^-1

#values kept constant in simulations
u = 1.66e-27 #1 atomic mass unit in kg

m_Mg = 24.31 #atomic mass magnesium
m_Si = 28.08 #atomic mass silisium
m_O = 16.00 #atomic mass oxygen
m_Fe = 55.84 #atomic mass iron
m_C = 12.01 #atomic mass carbon

M_ms = (m_Mg + m_Si + m_Fe + 4 * m_O) * 10**(-3) #molar mass silicate in kg per mol
M_mc = m_C * 10**(-3) #molar mass carbon in kg per mol

mA_S = (m_Mg + m_Si + m_Fe + 4 * m_O) / 7 * u #total mass, all constituents, in kg, silicate
mA_C = m_C * u #total mass in kg, carbon

R = 1 * au #initial radial position

"""scaled initial parameters"""
rhat0 = 1.0 #initial scaled radial position
thetahat0 = 0 #initial scaled angular position
vrhat0 = 0  #initial scaled radial velocity
omegahat0 = 1.0 #initial scaled angular velocity
betahat0 = 1.0 #initial scaled beta
mhat0 = 1.0 #initial scaled mass

init_polar_scaled = np.array([rhat0 , thetahat0 , vrhat0 , omegahat0]) #initial scaled polar coords
init_cart_scaled = polar_to_cartesian(init_polar_scaled) #initial scaled cart coords

sil_beta = "sil_radpr_prdrag_sun1au.dat" #Utilizing data sets provided by Li, A., 2026
car_beta = "ac_radpr_prdrag_sun1au.dat" #Utilizing data sets provided by Li, A., 2026

init_vals = {"large":{
            "r": 1.54079 * 10**(-6),
            "B": {"silicate": 0.1235 , "carbon": 0.2646} 
            },

            "medium":{
            "r": 0.17508 * 10**(-6),
            "B": {"silicate": 0.8560 , "carbon": 3.0589}
            },

            "small":{
            "r": 0.04259 * 10**(-6),
            "B":{"silicate": 0.2098 , "carbon": 1.6179}     
            },

            "biggest":{
                "r": 50.00000 * 10**(-6) ,
                "B": {"silicate": 0.0032 , "carbon": 0.0063}
            }}
            

material_files = {"silicate": sil_beta , "carbon": car_beta}

r_vals = np.linspace(1e-9 , 50e-6 , 200) #size range corresponding to that for realistic beta vals

yr = 60 * 60 * 24 * 365 #one year in s



"""converting .dat file to array with specified delimiter and column names"""
def dat_to_arr(file):
    """input: file (.dat)
     
       returns: dust_size , dust_betaval , dust_PRtime (array), containing dust sizes, beta values and Poynting-Robertson lifetime"""
    
    dust = pd.read_csv(file , sep = r'\s+' , header = None , names = ["Dust size (micron)" , "Beta value" , "PR drag time (years)"])
    dust_size = dust["Dust size (micron)"].to_numpy() 
    dust_betaval = dust["Beta value"].to_numpy()
    dust_PRtime = dust["PR drag time (years)"].to_numpy() 

    return dust_size , dust_betaval , dust_PRtime

def size_to_mass(r , material):
    if material == "silicate":
        rho = rho_s

    else: 
        rho = rho_c

    m = 4 / 3 * np.pi * rho * r**3

    return m

"""evaluating particle size against beta0"""
r_vals = np.linspace(0.00100 * 10**(-6) , 50 * 10**(-6) , 200)
r_betatest = np.linspace(0.00100 , 50 , 200)
m_range = size_to_mass(r_vals , "silicate") #masses corresponding to size range silicate
machine_eps = 1e-8

sil_size , sil_betaval , sil_PR = dat_to_arr(sil_beta) #fetching silicate size and beta values
sil_size = np.array(sil_size * 1e-6)
sil_betaval = np.array(sil_betaval)
sil_PR = np.array(sil_PR)

sil_mass = size_to_mass(sil_size , "silicate")
car_size , car_betaval , car_PR = dat_to_arr(car_beta) #fetching carbon size and beta values
car_size = np.array(car_size * 1e-6)
car_betaval = np.array(car_betaval)
car_PR = np.array(car_PR)
car_mass = size_to_mass(car_size , "carbon")

mask = (1 - car_betaval) > machine_eps

car_betaval_bound = car_betaval[mask]
car_size_bound = car_size[mask]
car_PR_bound = car_PR[mask]

car_mass_bound = size_to_mass(car_size_bound , "carbon")

material_files_bound = {"silicate": (sil_size , sil_betaval , sil_PR) , 
                        "carbon": (car_size_bound , car_betaval_bound , car_PR_bound)}

diffs = np.diff(sil_betaval[100:])
diffs_max = np.argmax(diffs)


"""t hat combinations used, dt timestep, t_tot total simulation time"""
#1000 orbits
dt5 = 3.16e3 #/ 6834848 #original T
t_tot5 = 1000
t5 = (dt5 , t_tot5)

#10000 orbits
dt6 = 3.16e3 #/ 6834848
t_tot6 = 10000
t6 = (dt6 , t_tot6)

#20000 orbits
dt7 = 3.16e3 #/ 6834848
t_tot7 = 20000
t7 = (dt7 , t_tot7)

dt8 = 3.16e3
t_tot8 = 12000
t8 = (dt8 , t_tot8)


if __name__ == "__main__":
    1
    print(4/3*np.pi*3500*(2.419656616e-13)**3)


