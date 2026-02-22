import numpy as np
from AstronomicalSilicate_modified import sputter

#values kept constant in simulations

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

mA_S = (m_Mg + m_Si + m_Fe + 4 * m_O) / 7 * u #total mass, all constituents, in kg, silicate
mA_C = m_C * u #total mass in kg, carbon


if __name__ == "__main__":
    a = sputtering_yield()
    print(a)


