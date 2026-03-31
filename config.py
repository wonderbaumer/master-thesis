import numpy as np
import pandas as pd

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
            }

        }

material_files = {"silicate": sil_beta , "carbon": car_beta}

r_vals = np.linspace(1e-9 , 50e-6 , 200) #size range corresponding to that for realistic beta vals

yr = 60 * 60 * 24 * 365 #one year in s

"""radiation pressure"""
rho = 2500 #kgm^-3
S_s = 1361 #solar constant, in Wm^-2
au = 149597871e3 #one astronomical unit, AU, in m
q_pr = 1 #radiation pressure coefficient, unitless
m_s = 1.98847e30  #mass of sun, in kg
c = 299792458 #speed of light, in ms^-1



"""converting .dat file to array with specified delimiter and column names"""
def dat_to_arr(file):
    """input: file (.dat)
     
       returns: dust_size , dust_betaval , dust_PRtime (array), containing dust sizes, beta values and Poynting-Robertson lifetime"""
    
    dust = pd.read_csv(file , sep = r'\s+' , header = None , names = ["Dust size (micron)" , "Beta value" , "PR drag time (years)"])
    dust_size = dust["Dust size (micron)"].to_numpy() 
    dust_betaval = dust["Beta value"].to_numpy()
    dust_PRtime = dust["PR drag time (years)"].to_numpy() 

    return dust_size , dust_betaval , dust_PRtime

"""evaluating particle size against beta0"""
r_vals = np.linspace(0.00100 * 10**(-6) , 50 * 10**(-6) , 200)
r_betatest = np.linspace(0.00100 , 50 , 200)
m_range = 4 / 3 * np.pi * rho * r_vals**3 #masses corresponding to size range

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


if __name__ == "__main__":
    sil_size , sil_betaval , sil_PRtime = dat_to_arr(sil_beta)
    fsw = sw_flux()
    Ytot = sputtering_yield("silicate" , "slow" , "all")
    


