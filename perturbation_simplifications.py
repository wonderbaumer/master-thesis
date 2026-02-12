import sympy as sp
from config import rhat0 , betahat0

"""defining the hatted variables"""
t = sp.Symbol("t") #time
r_0 = sp.Function("r_0")(t) #r0
r_10 = sp.Function("r_10")(t) #r10
r_01 = sp.Function("r_01")(t) #r01

v_r0 = sp.Function("v_r0")(t) #v0
v_r10 = sp.Function("v_r10")(t) #vr10
v_r01 = sp.Function("v_r01")(t) #vr01

vrdot_0 = sp.Function("vrdot_0")(t) #vdot0
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

beta_0 = sp.Symbol("beta_0") #beta0
beta_10 = sp.Symbol("beta_10") #beta10
beta_01 = sp.Symbol("beta_01") #beta01

B = sp.Symbol("B") #B, initial beta 
V = sp.Symbol("V") #V, initial velocity
c = sp.Symbol("c") #c, speed of light

epsilon = sp.Symbol("epsilon") #epsilon
delta = sp.Symbol("delta") #delta


r = r_0 + epsilon * r_10 + delta * r_01 #r perturbed expression
vr = v_r0 + epsilon * v_r10 + delta * v_r01 #v perturbed expression
vr_dot = vrdot_0 + epsilon * vrdot_10 + delta * vrdot_01 #vdot perturbed expression

beta = beta_0 + epsilon * beta_10 + delta * beta_01 #beta petrubed expression

theta = theta_0 + epsilon * theta_10 + delta * theta_01 #theta perturbed expression
omega = omega_0 + epsilon * omega_10 + delta * omega_01 #omega perturbed expression
omega_dot = omegadot_0 + epsilon * omegadot_10 + delta * omegadot_01 #omegadot perturbed expression

"""radial equation"""
rad_eq = sp.Eq((vr_dot - r * omega**2) * (1 - B) * r**2 , 
               (-1 + B * beta) * (1 - (2 * vr * delta))) #radial eq of motion

"""zeroth and first order expressions"""
req_lhs = sp.series(rad_eq.lhs , epsilon , 0 , 2).removeO() #removing O(epsilon^2) lhs
req_lhs = sp.series(req_lhs , delta , 0 , 2).removeO() #removing O(delta^2) lhs

req_rhs = sp.series(rad_eq.rhs , epsilon , 0 , 2).removeO() #removing O(epsilon^2) rhs
req_rhs = sp.series(req_rhs , delta , 0 , 2).removeO() #removing O(delta^2) rhs

rad_eq = (req_lhs - req_rhs).expand()
rad_eq_zeroth_order = rad_eq.coeff(epsilon , 0).coeff(delta , 0) #zeroth order total expression
rad_eq_10 = rad_eq.coeff(epsilon , 1).coeff(delta , 0) #10 expression
rad_eq_01 = rad_eq.coeff(epsilon , 0).coeff(delta , 1) #01 expression

vrdot0_sol = sp.solve(rad_eq_zeroth_order , vrdot_0) #zeroth order, solving for vdot0
vrdot10_sol = sp.solve(rad_eq_10 , vrdot_10) #first order, solving for vdot10
vrdot01_sol = sp.solve(rad_eq_01 , vrdot_01) #first order, solving for vdot01

"""angular equation"""
ang_eq = sp.Eq(r**2 * (1 - B) * (r * omega_dot + 2 * vr * omega) , (1 - B * beta) * r * omega * delta) #angular eq of motion

"""first order sols"""
angeq_lhs = sp.series(ang_eq.lhs , epsilon , 0 , 2).removeO() #removing O(epsilon^2) lhs
angeq_lhs = sp.series(angeq_lhs , delta , 0 , 2).removeO() #removing O(delta^2) lhs

angeq_rhs = sp.series(ang_eq.rhs , epsilon , 0 , 2).removeO() #removing O(epsilon^2) rhs
angeq_rhs = sp.series(angeq_rhs , delta , 0 , 2).removeO() #removing O(delta^2) rhs

angeq = (angeq_lhs - angeq_rhs).expand()
angeq_zeroth_order = angeq.coeff(epsilon , 0).coeff(delta , 0) #zeroth order total expression
angeq_10 = angeq.coeff(epsilon , 1).coeff(delta , 0) #10 expression
angeq_01 = angeq.coeff(epsilon , 0).coeff(delta , 1) #01 expression

omegadot_0 = sp.solve(angeq_zeroth_order , omegadot_0) #zeroth order, omegadot0
omegadot_10 = sp.solve(angeq_10 , omegadot_10) #first order, omegadot10
omegadot_01 = sp.solve(angeq_01 , omegadot_01) #first order, omegadot01

print(omegadot_10)