from scipy.constants import *
import numpy as np

yr = 60 * 60 * 24 * 365

"__radiation pressure__"
rho = 2500 #kgm^-3
S_s = 1361 #solar constant, in Wm^-2
au = 149597871e3 #one astronomical unit, AU, in m
q_pr = 1 #radiation pressure coefficient, unitless
m_s = 1.98847e30  #mass of sun, in kg

"__beta calcs___"
beta_const = S_s * q_pr * np.pi * au**2 / (G * m_s * c) #beta const
betahat_0 = 1 #initial betahat, when beta=beta_0

"__sputtering__"
N_sw = 3e6 #m^-3, density slow solar wind
v_sw = 800e3 #ms^-1, velocity slow solar wind
fsw = N_sw * v_sw #solar wind flux, m^-2s^-1

YH = 0.05 #sputtering yield hydrogen, slow sw
YHe = 0.02 #sputtering yield helium, slow sw
Ytot = YH + YHe #total sputtering yield, slow sw

u = 1.66e-27 #1 atomic mass unit in kg

m_Mg = 24.31 #atomic mass magnesium
m_Si = 28.08 #atomic mass silisium
m_O = 16.00 #atomic mass oxygen

mA = (m_Mg + m_Si + 4 * m_O) / 7 * u #total mass, all constituents, in kg




