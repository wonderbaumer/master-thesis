import sympy as sp

"""This code has perturbed the particle in one parameter, epsilon_0, and solves it algebraically"""

#Defining the variables
t = sp.Symbol("t") #time
dt = sp.Symbol("dt") #time differential

r_0 = sp.Function("r_0")(t) #r0
r_1 = sp.Function("r_1")(t) #r1

v_r0 = sp.Function("v_r0")(t) #vr0
v_r1 = sp.Function("v_r1")(t) #vr1

vrdot_0 = sp.Function("vrdot_0")(t) #vrdot0
vrdot_1 = sp.Function("vrdot_1")(t) #vrdot1

theta_0 = sp.Function("theta_0")(t) #theta0
theta_1 = sp.Function("theta_1")(t) #theta1

omega_0 = sp.Function("omega_0")(t) #omega0
omega_1 = sp.Function("omega_1")(t) #omega1

omegadot_0 = sp.Function("omegadot_0")(t) #omegadot0
omegadot_1 = sp.Function("omegadot_1")(t) #omegadot1

"""Beta and mass relations are derived by hand"""
beta_0 = sp.Function("beta0")(t) #beta0
beta_1 = sp.Function("beta1")(t) #beta1

dt_r0 = sp.Symbol("dt_r0") #dr0/dt
dt_r1 = sp.Symbol("dt_r1") #dr1/dt
dtt_r0 = sp.Symbol("dtt_r0") #d^2r0/dt^2
dtt_r1 = sp.Symbol("dtt_r1") #d^2r1/dt^2

dt_theta0 = sp.Symbol("dt_theta0") #dtheta0/dt
dt_theta1 = sp.Symbol("dt_theta1") #dtheta1/dt
dtt_theta0 = sp.Symbol("dtt_theta0") #d^2theta0/dt^2
dtt_theta1 = sp.Symbol("dtt_theta1") #d^2theta1/dt^2

B = sp.Symbol("B") #B, initial beta 
delta = sp.Symbol("delta") #delta, drag term

epsilon_0 = sp.Symbol("epsilon_0") #epsilon_0

r_exp = r_0 + epsilon_0 * r_1 #r perturbed expression
vr_exp = dt_r0 + epsilon_0 * dt_r1  #v perturbed expression
vrdot_exp = dtt_r0 + epsilon_0 * dtt_r1 + epsilon_0 * ?????  #vrdot perturbed expression

theta_exp = theta_0 + epsilon_0 * E_r * theta_1 #theta perturbed expression
omega_exp = dt_theta0 - 2 * epsilon_0 * theta_1 / r_0**3 * dt_r0 + epsilon_0 / r_0**2 * dt_theta1  #omega perturbed expression
omegadot_exp = dtt_theta0 - 2 * epsilon_0 * theta_1 / r_0**3 * dtt_r0  + omegadotterm #omegadot perturbed expression

beta_exp = beta_0 + epsilon_0 * E_r * beta_1 #beta perturbed expression

#Radial equation
rad_eq = sp.Eq((1 - B) * r_exp**2 * (vrdot_exp - r_exp * omega_exp**2) , 
               -(1 - beta_exp * B) - 2 * delta * beta_exp * B * vr_exp) #radial eq of motion

#up to second order expressions
req_lhs = sp.series(rad_eq.lhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) lhs

req_rhs = sp.series(rad_eq.rhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) rhs

rad_eq = (req_lhs - req_rhs).expand()
rad_eq_zeroth_order = rad_eq.coeff(epsilon_0 , 0) #zeroth order total expression
rad_eq_1 = rad_eq.coeff(epsilon_0 , 1) #1 order expression

#angular equation
ang_eq = sp.Eq(r_exp * (1 - B) * (r_exp * omegadot_exp + 2 * vr_exp * omega_exp) , -delta * B * beta_exp * omega_exp) #angular eq of motion

#second order sols
angeq_lhs = sp.series(ang_eq.lhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) lhs

angeq_rhs = sp.series(ang_eq.rhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) rhs

angeq = (angeq_lhs - angeq_rhs).expand()
angeq_zeroth_order = angeq.coeff(epsilon_0 , 0) #zeroth order total expression
angeq_1 = angeq.coeff(epsilon_0 , 1) #1 expression

print(rad_eq_zeroth_order)