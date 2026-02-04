import sympy as sp
from config import rhat0 , betahat0

"""defining the hatted variables"""
t = sp.Symbol("t") #time
r_0 = sp.Function("r_0")(t) #r0
r_1 = sp.Function("r_1")(t) #r1

v_r0 = sp.Function("v_r0")(t) #v0
v_r1 = sp.Function("v_r1")(t) #v1

vrdot_0 = sp.Function("vrdot_0")(t) #vdot0
vrdot_1 = sp.Function("vrdot_1")(t) #vdot1

theta_0 = sp.Function("theta_0")(t) #theta0
theta_1 = sp.Function("theta_1")(t) #theta1

omega_0 = sp.Function("omega_0")(t) #omega0
omega_1 = sp.Function("omega_1")(t) #omega1

omegadot_0 = sp.Function("omegadot_0")(t) #omegadot0
omegadot_1 = sp.Function("omegadot_1")(t) #omegadot1

beta_0 = sp.Symbol("beta_0") #beta0
B = sp.Symbol("B") #B, initial beta 
beta_1 = sp.Symbol("beta_1") #beta1
epsilon = sp.Symbol("epsilon") #epsilon
V = sp.Symbol("V") #V, initial velocity
c = sp.Symbol("c") #c, speed of light

r = r_0 + epsilon * r_1 #r perturbed expression
vr = v_r0 + epsilon * v_r1 #v perturbed expression
vr_dot = vrdot_0 + epsilon * vrdot_1 #vdot perturbed expression

beta = beta_0 + epsilon * beta_1 #beta petrubed expression

theta = theta_0 + epsilon * theta_1 #theta perturbed expression
omega = omega_0 + epsilon * omega_1 #omega perturbed expression
omega_dot = omegadot_0 + epsilon * omegadot_1 #omegadot perturbed expression

"""radial equation"""
rad_eq = sp.Eq((vr_dot - r * omega**2) * (1 - B) * r**2 , 
               (-1 + B * beta) * (1 - (2 * vr * V / c))) #radial eq of motion

"""zeroth and first order expressions"""
req_lhs = sp.series(rad_eq.lhs , epsilon , 0 , 2).removeO() #removing O(epsilon^2) lhs
req_rhs = sp.series(rad_eq.rhs , epsilon , 0 , 2).removeO() #removing O(epsilon^2) rhs

rad_eq_zeroth_order = (req_lhs - req_rhs).expand().coeff(epsilon , 0) #zeroth order total expression
rad_eq_first_order = (req_lhs - req_rhs).expand().coeff(epsilon , 1) #first order total expression

vrdot0_sol = sp.solve(rad_eq_zeroth_order , vrdot_0) #zeroth order, solving for vdot0
vrdot1_sol = sp.solve(rad_eq_first_order , vrdot_1) #first order, solving for vdot1

"""zeroth and first order solutions"""
vdot0_sol = vrdot0_sol[0].subs({
    beta_0 : betahat0 , 
    r_0 : rhat0
    })
vdot1_sol = vrdot1_sol[0].subs({
    beta_0 : betahat0 , 
    r_0 : rhat0
    })

"""angular equation"""
ang_eq = sp.Eq(r**2 * (1 - B) * (r * omega_dot + 2 * vr * omega) , (1 - B * beta) * r * omega * V / c) #angular eq of motion

"""first order sols"""
angeq_lhs = sp.series(ang_eq.lhs , epsilon , 0 , 2).removeO() #removing O(epsilon^2) lhs
angeq_rhs = sp.series(ang_eq.rhs , epsilon , 0 , 2).removeO() #removing O(epsilon^2) rhs

angeq_zeroth_order = (angeq_lhs - angeq_rhs).expand().coeff(epsilon , 0) #zeroth order total expression
angeq_first_order = (angeq_lhs - angeq_rhs).expand().coeff(epsilon , 1) #first order total expression

omegadot_0 = sp.solve(angeq_zeroth_order , omegadot_0) #zeroth order, omegadot0
omegadot_1 = sp.solve(angeq_first_order , omegadot_1) #zeroth order, omegadot1

omega_0 = sp.solve(angeq_zeroth_order , omega_0) #zeroth order, omega0
omega_1 = sp.solve(angeq_first_order , omega_1) #first order, omega1

print(omegadot_1)