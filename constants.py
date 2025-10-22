from scipy.constants import *
import numpy as np

rho = 2500 #kgm^-3
S_s = 1361 #solar constant, in Wm^-2
au = 149597871e3 #one astronomical unit, AU, in m
r0 = 1.0 * au #initial radial dist in units of AU
q_pr = 1 #radiation pressure coefficient, unitless
m_s = 1.98847e30  #mass of sun, in kg
yr = 60 * 60 * 24 * 365

beta_const = S_s * q_pr * np.pi * au**2 / (G * m_s * c) #beta const
beta_0 = 0.45931933916320633 #initial mass and radius in cst table

period_cst = np.sqrt(r0**3 / (G * m_s))

N_sw = 3e6 #m^-3, density slow solar wind
v_sw = 800e3 #ms^-1, velocity slow solar wind
fsw = N_sw * v_sw #solar wind flux, m^-2s^-1

YH = 0.05 #sputtering yield hydrogen
YHe = 0.02 #sputtering yield helium
Ytot = YH + YHe #total sputtering yield

u = 1.66e-24 #1 atomic mass unit in g

m_Mg = 24.31 #atomic mass
m_Si = 28.08 
m_O = 16.00

mA = (m_Mg + m_Si + 4 * m_O) / 7 * u

mass_cst = fsw * Ytot * mA * np.pi * (3 / (4 * np.pi * rho))**(2 / 3) #mass constant
