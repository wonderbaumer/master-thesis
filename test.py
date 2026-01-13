import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import G

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
rho = 2500 #kgm^-3
m_par = 1.30899694e-15 #mass particle in kg

B = 0.45931933916320633 #initial mass and radius from initial pressure radiation/gravity
au = 149597871e3 #one astronomical unit, AU, in m
m_s = 1.98847e30  #mass of sun, in kg
R = au #initial radial position
T = round(np.sqrt(R**3 / (G * m_s * (1 - B)))) #initial period, scaled formula


#500 orbits
dt4 = 3.16e3 #old e3 new e2
t_tot4 = 500 * T
t4 = (dt4 , t_tot4)

yr = 60 * 60 * 24 * 365 #one year in s

#temp500.
dt4 = 3.16e3
t_tot4 = 500 * 1.36 * yr
t4 = (dt4 , t_tot4)


rk = np.load("C:/Users/cecil/Documents/Project-paper/Files/rk45_t4_massloss_yrs.npz")
x_r = rk["x"]
y_r = rk["y"]
vx_r = rk["vx"]
vy_r = rk["vy"]
m_r = rk["m"]
b_r = rk["b"]

num_angle = np.atan2(y_r , x_r)
num_angle = np.unwrap(num_angle)

dt , t_tot = t4
t_hat = np.arange(0 , t_tot , dt) / T

def theta(t):
    f = B / (3 * (1 - B))
    theta0 = t
    theta1 = -f * t**2 - 2 * f * np.cos(t)

    epsilon = fsw * Ytot * mA * np.pi * (3 / (4 * np.pi * rho))**(2 / 3) * m_par**(-1 / 3) * T

    theta = theta0 + epsilon * theta1

    return theta

pert_angle = theta(t_hat)

print(num_angle , pert_angle)