import sympy as sp
from sympy import pprint , abc , dsolve , Derivative , integrate
from sympy.abc import a , b , c , g , e , d , f , h , i , j , k , l , m , n , o , p , q , r
from sympy.solvers.ode.systems import dsolve_system
from config import rhat0 , betahat0
from sympy.matrices.expressions import MatMul


#defining the hatted variables
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

beta_0 = sp.Function("beta0")(t)
beta_10 = sp.Function("beta10")(t)
beta_01 = sp.Function("beta01")(t)

dt_r0 = sp.Symbol("dt_r0")
dt_r10 = sp.Symbol("dt_r10")
dt_r01 = sp.Symbol("dt_r01")
dtt_r0 = sp.Symbol("dtt_r0")
dtt_r10 = sp.Symbol("dtt_r10")
dtt_r01 = sp.Symbol("dtt_r01")

dt_theta0 = sp.Symbol("dt_theta0")
dt_theta10 = sp.Symbol("dt_theta10")
dt_theta01 = sp.Symbol("dt_theta01")
dtt_theta0 = sp.Symbol("dtt_theta0")
dtt_theta10 = sp.Symbol("dtt_theta10")
dtt_theta01 = sp.Symbol("dtt_theta01")

B = sp.Symbol("B") #B, initial beta 
delta = sp.Symbol("delta")

epsilon_0 = sp.Symbol("epsilon_0") #epsilon_0
E_r = sp.Function("E")(r_0) #radial dependent epsilon

vrdotterm = - 2 * r_0**(-3) * epsilon_0 * dtt_r0 * r_10 + 6 * r_0**(-4) * epsilon_0 * dt_r0**2 * r_10 + delta * dtt_r01 - 2 * r_0**(-3) * epsilon_0 * dt_r0 * dt_r10
omegaterm = - 2 * r_0**(-3) * epsilon_0 * dt_r0 * theta_10 + delta * dt_theta01
omegadotterm = - 2 * epsilon_0 * r_0**(-3) * dt_r0 * dt_theta10 - 2 * r_0**(-3) * epsilon_0 * dt_r0 * dt_theta10 - 2 * r_0**(-3) * epsilon_0 * dtt_r0 * theta_10 + 6 * r_0**(-4) * epsilon_0 * dt_r0**2 * theta_10 + delta * dtt_theta01 

r_exp = r_0 + epsilon_0 * E_r * r_10 + delta * r_01 #r perturbed expression
vr_exp = dt_r0 + dt_r10 * r_0**(-2) * epsilon_0 - 2 * r_0**(-3) * epsilon_0 * dt_r0 + delta * dt_r01  #vr perturbed expression
vrdot_exp = dtt_r0 + r_0**(-2) * epsilon_0 * dtt_r10 - 2 * r_0**(3) * epsilon_0 * dt_r10 * dt_r0 + vrdotterm  #vrdot perturbed expression

theta_exp = theta_0 + epsilon_0 * E_r * theta_10 + delta * theta_01 #theta perturbed expression
omega_exp = dt_theta0 + epsilon_0 * r_0**(-2) * dt_theta10 + omegaterm  #omega perturbed expression
omegadot_exp = dtt_theta0 + epsilon_0 * r_0**(-2) * dtt_theta10 + omegadotterm #omegadot perturbed expression

beta_exp = beta_0 + epsilon_0 * E_r * beta_10 + delta * beta_01 #beta perturbed expression

#Radial equation
rad_eq = sp.Eq((1 - B) * r_exp**2 * (vrdot_exp - r_exp * omega_exp**2) , 
               -(1 - beta_exp * B) - 2 * delta * beta_exp * B * vr_exp) #radial eq of motion

#up to second order expressions
req_lhs = sp.series(rad_eq.lhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) lhs
req_lhs = sp.series(req_lhs , delta , 0 , 2).removeO()

req_rhs = sp.series(rad_eq.rhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) rhs
req_rhs = sp.series(req_rhs , delta , 0 , 2).removeO()

rad_eq = (req_lhs - req_rhs).expand()
rad_eq_zeroth_order = rad_eq.coeff(epsilon_0 , 0) #zeroth order total expression
rad_eq_zeroth_order = rad_eq_zeroth_order.coeff(delta , 0)

rad_eq_10 = rad_eq.coeff(epsilon_0 , 1).coeff(delta , 0) #10 order expression
rad_eq_01 = rad_eq.coeff(epsilon_0 , 0).coeff(delta , 1)

#angular equation
ang_eq = sp.Eq(r_exp * (1 - B) * (r_exp * omegadot_exp + 2 * vr_exp * omega_exp) , -delta * B * beta_exp * omega_exp) #angular eq of motion

#second order sols
angeq_lhs = sp.series(ang_eq.lhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) lhs
angeq_lhs = sp.series(angeq_lhs , delta , 0 , 2).removeO()

angeq_rhs = sp.series(ang_eq.rhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) rhs
angeq_rhs = sp.series(angeq_rhs , delta , 0 , 2).removeO()

angeq = (angeq_lhs - angeq_rhs).expand()
angeq_zeroth_order = angeq.coeff(epsilon_0 , 0) #zeroth order total expression
angeq_zeroth_order = angeq_zeroth_order.coeff(delta , 0)

angeq_10 = angeq.coeff(epsilon_0 , 1).coeff(delta , 0) #1 expression
angeq_01 = angeq.coeff(epsilon_0 , 0).coeff(delta , 1)

print(angeq_01)