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
r_1 = sp.Function("r_1")(t) #r1

v_r0 = sp.Function("v_r0")(t) #v0
v_r1 = sp.Function("v_r1")(t) #vr1

vrdot_0 = sp.Function("vrdot_0")(t) #vdot0
vrdot_1 = sp.Function("vrdot_1")(t) #vrdot1

theta_0 = sp.Function("theta_0")(t) #theta0
theta_1 = sp.Function("theta_1")(t) #theta1

omega_0 = sp.Function("omega_0")(t) #omega0
omega_1 = sp.Function("omega_1")(t) #omega1

omegadot_0 = sp.Function("omegadot_0")(t) #omegadot0
omegadot_1 = sp.Function("omegadot_1")(t) #omegadot1

beta_0 = sp.Function("beta0")(t)
beta_1 = sp.Function("beta1")(t)

dt_r0 = sp.Symbol("dt_r0")
dt_r1 = sp.Symbol("dt_r1")
dtt_r0 = sp.Symbol("dtt_r0")
dtt_r1 = sp.Symbol("dtt_r1")

dt_theta0 = sp.Symbol("dt_theta0")
dt_theta1 = sp.Symbol("dt_theta1")
dtt_theta0 = sp.Symbol("dtt_theta0")
dtt_theta1 = sp.Symbol("dtt_theta1")

B = sp.Symbol("B") #B, initial beta 
V = sp.Symbol("V")
c = sp.Symbol("c")

epsilon_0 = sp.Symbol("epsilon_0") #epsilon_0
E_r = sp.Function("E")(r_0) #radial dependent epsilon
vrdotterm = 6 * epsilon_0 / r_0**4 * r_1 * dt_r0**2 + epsilon_0 / r_0**2 * dtt_r1 - 2 * epsilon_0 / r_0**3 * r_1 * dt_r0
omegadotterm = - 2 * epsilon_0 / r_0**3 * dt_r0 * dt_theta1 - 6 * epsilon_0 * theta_1 / r_0**4 * dt_r0**2 + epsilon_0 / r_0**2 * dtt_theta1 - 2 * epsilon_0 / r_0**3 * theta_1 * dt_r0 

r_exp = r_0 + epsilon_0 * E_r * r_1 #r perturbed expression
vr_exp = dt_r0 - 2 * epsilon_0 * r_1 / r_0**3 * dt_r0 + epsilon_0 / r_0**2 * dt_r1 #v perturbed expression
vrdot_exp = dtt_r0 * (1 - 2 * epsilon_0 * r_1 / r_0**3) - 2 * epsilon_0 / r_0**3 * dt_r0 * dt_r1 + vrdotterm  #vrdot perturbed expression

theta_exp = theta_0 + epsilon_0 * E_r * theta_1 #theta perturbed expression
omega_exp = dt_theta0 - 2 * epsilon_0 * theta_1 / r_0**3 * dt_r0 + epsilon_0 / r_0**2 * dt_theta1  #omega perturbed expression
omegadot_exp = dtt_theta0 - 2 * epsilon_0 * theta_1 / r_0**3 * dtt_r0  + omegadotterm #omegadot perturbed expression

beta_exp = beta_0 + epsilon_0 * E_r * beta_1 #beta perturbed expression

#Radial equation
rad_eq = sp.Eq((1 - B) * r_exp**2 * (vrdot_exp - r_exp * omega_exp**2) , 
               -(1 - beta_exp * B) - 2 * V / c * beta_exp * B * vr_exp) #radial eq of motion

#up to second order expressions
req_lhs = sp.series(rad_eq.lhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) lhs

req_rhs = sp.series(rad_eq.rhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) rhs

rad_eq = (req_lhs - req_rhs).expand()
rad_eq_zeroth_order = rad_eq.coeff(epsilon_0 , 0) #zeroth order total expression
rad_eq_1 = rad_eq.coeff(epsilon_0 , 1) #1 order expression

#angular equation
ang_eq = sp.Eq(r_exp * (1 - B) * (r_exp * omegadot_exp + 2 * vr_exp * omega_exp) , -V / c * B * beta_exp * omega_exp) #angular eq of motion

#second order sols
angeq_lhs = sp.series(ang_eq.lhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) lhs

angeq_rhs = sp.series(ang_eq.rhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) rhs

angeq = (angeq_lhs - angeq_rhs).expand()
angeq_zeroth_order = angeq.coeff(epsilon_0 , 0) #zeroth order total expression
angeq_1 = angeq.coeff(epsilon_0 , 1) #1 expression

print(rad_eq_zeroth_order)