import numpy as np
from AstronomicalSilicate_modified import sputter
import pandas as pd
from scipy.constants import N_A

#values kept constant in simulations

sil_beta = "sil_radpr_prdrag_sun1au.dat" #Utilizing data sets provided by Li, A., 2026
car_beta = "ac_radpr_prdrag_sun1au.dat" #Utilizing data sets provided by Li, A., 2026

r_vals = np.linspace(1e-9 , 50e-6 , 200) #size range corresponding to that for realistic beta vals

yr = 60 * 60 * 24 * 365 #one year in s

"""radiation pressure"""
rho = 2500 #kgm^-3
S_s = 1361 #solar constant, in Wm^-2
au = 149597871e3 #one astronomical unit, AU, in m
q_pr = 1 #radiation pressure coefficient, unitless
m_s = 1.98847e30  #mass of sun, in kg
c = 299792458 #speed of light, in ms^-1

#Solar wind flux
def sw_flux(sw = "slow"):
    """input: sw (string), solar wind conditions, default slow, options fast, slow, CME
        
        returns: fsw (array), solar wind flux"""
    
    if sw == "fast":
        N_sw = 3e6 #m^-3 density fast solar wind
        v_sw = 8e5 #ms^-1, velocity fast solar wind
        fsw = N_sw * v_sw #solar wind flux, m^-2s^-1
    
    elif sw == "slow":
        N_sw = 8e6 #m^-3, density slow solar wind
        v_sw = 3e5 #ms^-1, velocity slow solar wind
        fsw = N_sw * v_sw #solar wind flux, m^-2s^-1
    
    elif sw == "CME":
         N_sw = 7e7 #m^-3, density slow solar wind
         v_sw = 5e5 #ms^-1, velocity slow solar wind
         fsw = N_sw * v_sw #solar wind flux, m^-2s^-1
    
    return fsw

#Total sputtering yield
def sputtering_yield(material = "silicate" , sw = "slow" , species = "all"):
    """input: sw (string), solar wind conditions, default slow, options fast, slow, CME
              species (string), solar wind elements, default all, else one of H, He, C, O, N, Fe, Ne, Mg, Si, S
        
        returns: Ytot (array), total sputtering yield"""
    
    if species == "all":
        H = np.sum(sputter[material][sw]["H"]) #total sputtering yield hydrogen
        He = np.sum(sputter[material][sw]["He"]) #total sputtering yield helium
        C = np.sum(sputter[material][sw]["C"]) #total sputtering yield carbon
        O = np.sum(sputter[material][sw]["O"]) #total sputtering yield oxygen
        N = np.sum(sputter[material][sw]["N"]) #total sputtering yield nitrogen
        Fe = np.sum(sputter[material][sw]["Fe"]) #total sputtering yield iron
        Ne = np.sum(sputter[material][sw]["Ne"]) #total sputtering yield neon
        Mg = np.sum(sputter[material][sw]["Mg"]) #total sputtering yield magnesium
        Si = np.sum(sputter[material][sw]["Si"]) #total sputtering yield silicon
        S = np.sum(sputter[material][sw]["S"]) #total sputtering yield sulfur

        Ytot = H + He + C + O + N + Fe + Ne + Mg + Si + S #total sputtering yield 
    
    else:
        spec = np.sum(sputter[material][sw][species]) #total sputtering yield specified species

        Ytot = spec #total sputtering yield

    return Ytot

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

"""sputtering lifetime calcs"""
def sputtering_lifetime(r0 , fsw , Ytot , M_m):
    """input: r0 (float), initial dust radius
              fsw (float), solar wind flux
              Ytot (float), total sputtering yield
              M_m (float), molar mass of the material

       returns: t_sp (float), sputtering lifetime"""

    t_sp = (4 * r0 * N_A * rho) / (fsw * Ytot * M_m) #sputtering lifetime
    t_sp = t_sp / yr

    return t_sp

"""converting .dat file to array with specified delimiter and column names"""
def dat_to_arr(file):
    """input: file (.dat)
     
       returns: dust_size , dust_betaval , dust_PRtime (array), containing dust sizes, beta values and Poynting-Robertson lifetime"""
    
    dust = pd.read_csv(file , sep = r'\s+' , header = None , names = ["Dust size (micron)" , "Beta value" , "PR drag time (years)"])
    dust_size = dust["Dust size (micron)"].to_numpy() 
    dust_betaval = dust["Beta value"].to_numpy()
    dust_PRtime = dust["PR drag time (years)"].to_numpy() 

    return dust_size , dust_betaval , dust_PRtime


if __name__ == "__main__":
    sil_size , sil_betaval , sil_PRtime = dat_to_arr(sil_beta)
    fsw = sw_flux()
    Ytot = sputtering_yield("silicate" , "slow" , "all")
    


