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

"""sputtering"""
N_sw = 8e6 #m^-3, density slow solar wind
v_sw = 3e5 #ms^-1, velocity slow solar wind
fsw = N_sw * v_sw #solar wind flux, m^-2s^-1

#Slow solar wind, sputtering yields
def sputtering_yield(sw , species):
    if species == "all":

        H = np.sum(sputter[sw]["H"])
        He = np.sum(sputter[sw]["He"])
        C = np.sum(sputter[sw]["C"])
        O = np.sum(sputter[sw]["O"])
        N = np.sum(sputter[sw]["N"])
        Fe = np.sum(sputter[sw]["Fe"])
        Ne = np.sum(sputter[sw]["Ne"])
        Mg = np.sum(sputter[sw]["Mg"])
        Si = np.sum(sputter[sw]["Si"])
        S = np.sum(sputter[sw]["S"])

        Ytot = H + He + C + O + N + Fe + Ne + Mg + Si + S
    
    else:
        spec = np.sum(sputter[sw][species])
        Ytot = spec

    return Ytot

u = 1.66e-27 #1 atomic mass unit in kg

m_Mg = 24.31 #atomic mass magnesium
m_Si = 28.08 #atomic mass silisium
m_O = 16.00 #atomic mass oxygen

mA = (m_Mg + m_Si + 4 * m_O) / 7 * u #total mass, all constituents, in kg


if __name__ == "__main__":
    a = sputtering_yield("slow" , "all")
    print(a)


