import sympy as sp

"""Defining the hatted variables"""
t = sp.Symbol("t") #time
dt = sp.Symbol("dt") #time differential

r_0 = sp.Function("r_0")(t) #r0
r_10 = sp.Function("r_10")(t) #r10
r_01 = sp.Function("r_01")(t) #r01

v_r0 = sp.Function("v_r0")(t) #vr0
v_r10 = sp.Function("v_r10")(t) #vr10
v_r01 = sp.Function("v_r01")(t) #vr01

vrdot_0 = sp.Function("vrdot_0")(t) #vrdot0
vrdot_10 = sp.Function("vrdot_10")(t) #vrdot10
vrdot_01 = sp.Function("vrdot_01")(t) #vrdot01

theta_0 = sp.Function("theta_0")(t) #theta0
theta_10 = sp.Function("theta_10")(t) #theta10
theta_01 = sp.Function("theta_01")(t) #theta01

omega_0 = sp.Function("omega_0")(t) #omega0
omega_10 = sp.Function("omega_10")(t) #omega10
omega_01 = sp.Function("omega_01")(t) #omega01

omegadot_0 = sp.Function("omegadot_0")(t) #omegadot0
omegadot_10 = sp.Function("omegadot_10")(t) #omegadot10
omegadot_01 = sp.Function("omegadot_01")(t) #omegadot01

beta_0 = sp.Function("beta0")(t) #beta0
beta_10 = sp.Function("beta10")(t) #beta10
beta_01 = sp.Function("beta01")(t) #beta01

dt_r0 = sp.Symbol("dt_r0") #dr0/dt
dt_r10 = sp.Symbol("dt_r10") #dr10/dt
dt_r01 = sp.Symbol("dt_r01") #dr01/dt
dtt_r0 = sp.Symbol("dtt_r0") #d^2r0/dt^2
dtt_r10 = sp.Symbol("dtt_r10") #d^2r10/dt^2
dtt_r01 = sp.Symbol("dtt_r01") #d^2r01/dt^2

dt_theta0 = sp.Symbol("dt_theta0") #dtheta0/dt
dt_theta10 = sp.Symbol("dt_theta10") #dtheta10/dt
dt_theta01 = sp.Symbol("dt_theta01") #dtheta01/dt
dtt_theta0 = sp.Symbol("dtt_theta0") #d^2theta0/dt^2
dtt_theta10 = sp.Symbol("dtt_theta10") #d^2theta10/dt^2
dtt_theta01 = sp.Symbol("dtt_theta01") #d^2theta01/dt^2

B = sp.Symbol("B") #B, initial beta 
delta = sp.Symbol("delta") #delta, drag term
epsilon_0 = sp.Symbol("epsilon_0") #epsilon_0, mass loss rate   

r_exp = r_0 + epsilon_0 * r_10 + delta * r_01 #r perturbed expression
vr_exp = dt_r0 + dt_r10 * epsilon_0 + delta * dt_r01  #vr perturbed expression
vrdot_exp = dtt_r0 + epsilon_0 * dtt_r10 + delta * dtt_r01  #vrdot perturbed expression

theta_exp = theta_0 + epsilon_0 * theta_10 + delta * theta_01 #theta perturbed expression
omega_exp = dt_theta0 + epsilon_0 * dt_theta10 + delta * dt_theta01  #omega perturbed expression
omegadot_exp = dtt_theta0 + epsilon_0 * dtt_theta10 + delta * dtt_theta01  #omegadot perturbed expression

beta_exp = beta_0 + epsilon_0 * beta_10 + delta * beta_01 #beta perturbed expression

"""Radial equation"""
rad_eq = sp.Eq((1 - B) * r_exp**2 * (vrdot_exp - r_exp * omega_exp**2) , 
               -(1 - beta_exp * B) - 2 * delta * beta_exp * B * vr_exp) 

"""Up to epsilon^1 and delta^2 solutions"""
req_lhs = sp.series(rad_eq.lhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) lhs
req_lhs = sp.series(req_lhs , delta , 0 , 2).removeO() #removing O(delta^2) lhs

req_rhs = sp.series(rad_eq.rhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) rhs
req_rhs = sp.series(req_rhs , delta , 0 , 2).removeO() #removing O(delta^2) rhs

"""Separating into zeroth and first order in epsilon and delta solution"""
rad_eq = (req_lhs - req_rhs).expand()
rad_eq_zeroth_order = rad_eq.coeff(epsilon_0 , 0) #zeroth order in epsilon
rad_eq_zeroth_order = rad_eq_zeroth_order.coeff(delta , 0) #zeroth order in delta and epsilon

rad_eq_10 = rad_eq.coeff(epsilon_0 , 1).coeff(delta , 0) #first order in epsilon
rad_eq_01 = rad_eq.coeff(epsilon_0 , 0).coeff(delta , 1) #first order in delta

"""Angular equation"""
ang_eq = sp.Eq(r_exp * (1 - B) * (r_exp * omegadot_exp + 2 * vr_exp * omega_exp) 
               , -delta * B * beta_exp * omega_exp)

"""Up to epsilon^1 and delta^2 solutions"""
angeq_lhs = sp.series(ang_eq.lhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) lhs
angeq_lhs = sp.series(angeq_lhs , delta , 0 , 2).removeO() #removing O(delta^2) lhs

"""Separating into zeroth and first order in epsilon and delta solution"""
angeq_rhs = sp.series(ang_eq.rhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) rhs
angeq_rhs = sp.series(angeq_rhs , delta , 0 , 2).removeO() #removing O(delta^2) rhs

angeq = (angeq_lhs - angeq_rhs).expand()
angeq_zeroth_order = angeq.coeff(epsilon_0 , 0) #zeroth order in epsilon
angeq_zeroth_order = angeq_zeroth_order.coeff(delta , 0) #zeroth order in delta and epsilon

angeq_10 = angeq.coeff(epsilon_0 , 1).coeff(delta , 0) #first order in epsilon
angeq_01 = angeq.coeff(epsilon_0 , 0).coeff(delta , 1) #first order in delta

