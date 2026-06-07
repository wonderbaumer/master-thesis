import numpy as np
import pandas as pd
from polar_to_cart import polar_to_cartesian
from scipy.constants import c , G

#Radiation pressure
rho_s = 3500 #silicate bulk density, in kgm^-3
rho_c = 1800 #carbon bulk density, in kgm^-3
S_s = 1361 #solar constant, in Wm^-2
au = 149597871e3 #one astronomical unit, AU, in m
q_pr = 1 #radiation pressure coefficient, unitless
m_s = 1.98847e30  #mass of sun, in kg

#Values kept constant in simulations
u = 1.66e-27 #1 atomic mass unit, in kg

m_Mg = 24.31 #atomic mass magnesium, in u
m_Si = 28.08 #atomic mass silisium, in u
m_O = 16.00 #atomic mass oxygen, in u
m_Fe = 55.84 #atomic mass iron, in u
m_C = 12.01 #atomic mass carbon, in u

M_ms = (m_Mg + m_Si + m_Fe + 4 * m_O) * 10**(-3) #molar mass silicate, in kg per mol
M_mc = m_C * 10**(-3) #molar mass carbon, in kg per mol

mA_S = (m_Mg + m_Si + m_Fe + 4 * m_O) / 7 * u #mean mass, all constituents, in kg, silicate 
mA_C = m_C * u #mean mass, all constituens, in kg, carbon

R = 1 * au #initial radial position, in AU

#Scaled initial parameters
rhat0 = 1.0 #initial scaled radial position, unitless
thetahat0 = 0 #initial scaled angular position, unitless
vrhat0 = 0  #initial scaled radial velocity, unitless
omegahat0 = 1.0 #initial scaled angular velocity, unitless
betahat0 = 1.0 #initial scaled beta, unitless
mhat0 = 1.0 #initial scaled mass, unitless

init_polar_scaled = np.array([rhat0 , thetahat0 , vrhat0 , omegahat0]) #initial scaled polar coords, unitless
init_cart_scaled = polar_to_cartesian(init_polar_scaled) #initial scaled cart coords, unitless

sil_beta = "sil_radpr_prdrag_sun1au.dat" #Utilizing data sets provided by Li, A., 2026
car_beta = "ac_radpr_prdrag_sun1au.dat" #Utilizing data sets provided by Li, A., 2026

#Chosen initial values and labels
init_vals = {"A":{
            "r": 1.54079 * 10**(-6) ,
            "B": {"silicate": 0.1235 , "carbon": 0.2646} 
            },

            "C":{
            "r": 0.17508 * 10**(-6) ,
            "B": {"silicate": 0.8560 , "carbon": 3.0589}
            },

            "F":{
            "r": 0.04259 * 10**(-6) ,
            "B":{"silicate": 0.2098 , "carbon": 1.6179}     
            },

            "G": {
                "r": 0.01220 * 10**(-6) ,
                "B": {"silicate": 0.0832 , "carbon": 0.9758}
            } ,

            "H" : {
                "r" : 0.00708 * 10**(-6) ,
                "B" : {"silicate" : 0.0799 , "carbon" : 0.9505}
            } ,

            "D" : {
                "r" : 0.10165 * 10**(-6) ,
                "B" : {"silicate" : 0.6811 , "carbon" : 3.2628}
            } , 

            "E" : {
                "r" : 0.07745 * 10**(-6) ,
                "B" : {"silicate" : 0.5123 , "carbon" : 2.8217}
            }}
         
material_files = {"silicate": sil_beta , "carbon": car_beta} #Mapping labels to correct beta values

yr = 60 * 60 * 24 * 365 #one year, in s

"""converting .dat file to array with specified delimiter and column names"""
def dat_to_arr(file):
    """input: file (.dat)
     
       returns: dust_size , dust_betaval , dust_PRtime (array), containing dust sizes, beta values and Poynting-Robertson lifetime"""
    
    dust = pd.read_csv(file , sep = r'\s+' , header = None , names = ["Dust size (micron)" , "Beta value" , "PR drag time (years)"]) #unpacking file
    dust_size = dust["Dust size (micron)"].to_numpy() #column to array
    dust_betaval = dust["Beta value"].to_numpy() #column to array
    dust_PRtime = dust["PR drag time (years)"].to_numpy() #column to array

    return dust_size , dust_betaval , dust_PRtime

"""Converting from radius in m to mass in kg"""
def size_to_mass(r , material):
    """input: r (float), particle radius in m
              material (string), "silicate" or "carbon", the particle considered
        
       output: m (float), particle mass in kg
    """

    if material == "silicate":
        rho = rho_s #silicate bulk density

    else: 
        rho = rho_c #carbon bulk density

    m = 4 / 3 * np.pi * rho * r**3 #mass calcs based on ideal sphere

    return m

#Fetching particle sizes and beta values corresponding to bound orbit
tol = 1e-8

##Silicate
sil_size , sil_betaval , sil_PR = dat_to_arr(sil_beta) #fetching size, beta values and PR lifetime
sil_size = np.array(sil_size * 1e-6) #size from microns to metre
sil_betaval = np.array(sil_betaval) #beta values
sil_PR = np.array(sil_PR) #PR lifetime
sil_mass = size_to_mass(sil_size , "silicate") #size to mass

##Carbon
car_size , car_betaval , car_PR = dat_to_arr(car_beta) #fetching size, beta values and PR lifetime
car_size = np.array(car_size * 1e-6) #size from microns to metre
car_betaval = np.array(car_betaval) #beta values
car_PR = np.array(car_PR) #PR lifetime
car_mass = size_to_mass(car_size , "carbon") #size to mass

mask = (1 - car_betaval) > tol #masking beta values 1 or larger

car_betaval_bound = car_betaval[mask] #bound beta values
car_size_bound = car_size[mask] #bound sizes
car_PR_bound = car_PR[mask] #bound PR lifetime

car_mass_bound = size_to_mass(car_size_bound , "carbon") #bound masses

material_files_bound = {"silicate": (sil_size , sil_betaval , sil_PR) , 
                        "carbon": (car_size_bound , car_betaval_bound , car_PR_bound)} #bound materials mapping


#t hat combinations used, dt timestep, t_tot total simulation time
#1000 orbits
dt5 = 3.16e-3 
t_tot5 = 1000
t5 = (dt5 , t_tot5)

#10000 orbits
dt6 = 3.16e-3 
t_tot6 = 10000
t6 = (dt6 , t_tot6)

#20000 orbits
dt7 = 3.16e-3 
t_tot7 = 20000
t7 = (dt7 , t_tot7)

#30000 orbits
dt8 = 3.16e-3
t_tot8 = 30000
t8 = (dt8 , t_tot8) #la stå

#50000
dt10 = 3.16e-3
t_tot10 = 50000
t10 = (dt10 , t_tot10)

#100000
dt9 = 3.16e-3
t_tot9 = 100000
t9 = (dt9 , t_tot9)

def pr_lifetime(betaval , init_dist = 1.0):
        RAU = init_dist

        tau = 400 * RAU**2 / betaval

        return tau

tau_sil = pr_lifetime(sil_betaval)
tau_car = pr_lifetime(car_betaval)

if __name__ == "__main__":
    print(len(sil_size) , len(car_size_bound))
    
    

