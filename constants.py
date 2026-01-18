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

YH = 0.05 #sputtering yield hydrogen, slow sw
YHe = 0.02 #sputtering yield helium, slow sw
Ytot = YH + YHe #total sputtering yield, slow sw

u = 1.66e-27 #1 atomic mass unit in kg

m_Mg = 24.31 #atomic mass magnesium
m_Si = 28.08 #atomic mass silisium
m_O = 16.00 #atomic mass oxygen

mA = (m_Mg + m_Si + 4 * m_O) / 7 * u #total mass, all constituents, in kg




