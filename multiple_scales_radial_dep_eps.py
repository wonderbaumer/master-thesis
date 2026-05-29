import sympy as sp
from sympy import symbols , abc , dsolve , Derivative , integrate
from sympy.abc import a , b , c , g , e , d , f , h , i , j , k , l , m , n , o , p , q , r
from sympy.solvers.ode.systems import dsolve_system
from config import rhat0 , betahat0
from sympy.matrices.expressions import MatMul

#defining the hatted variables
t0 = sp.Symbol("t0") #fast time
t1 = sp.Symbol("t1") #slow time, symbol only
dt0 = sp.Symbol("dt0") #partial differential, fast time
dt1 = sp.Symbol("dt1") #partial differential, slow time

r_0 = sp.Function("r_0")(t1) #r0
r_1 = sp.Function("r_1")(t0 , t1) #r1

v_r0 = sp.Function("v_r0")(t1) #v0
v_r1 = sp.Function("v_r1")(t0 , t1) #vr1

vrdot_0 = sp.Function("vrdot_0")(t1) #vdot0
vrdot_1 = sp.Function("vrdot_1")(t0 , t1) #vrdot1

theta_0 = sp.Function("theta_0")(t0 , t1) #theta0
theta_1 = sp.Function("theta_1")(t0 , t1) #theta1

omega_0 = sp.Function("omega_0")(t1) #omega0
omega_1 = sp.Function("omega_1")(t0 , t1) #omega1

omegadot_0 = sp.Function("omegadot_0")(t1) #omegadot0
omegadot_1 = sp.Function("omegadot_1")(t0 , t1) #omegadot1

m_0 = sp.Function("m_0")(t1) #m0
m_1 = sp.Function("m_1")(t0 , t1) #m1

dt0_r0 = sp.Symbol("dt0_r0")
dt0_r1 = sp.Symbol("dt0_r1")
dt1_r0 = sp.Symbol("dt1_r0")
dt0t0_r0 = sp.Symbol("dt0t0_r0")
dt0t0_r1 = sp.Symbol("dt0t0_r1")
dt0t1_r0 = sp.Symbol("dt0t1_r0")
dt1t0_r0 = sp.Symbol("dt1t0_r0")

dt0_theta0 = sp.Symbol("dt0_theta0")
dt1_theta0 = sp.Symbol("dt1_theta0")
dt0_theta1 = sp.Symbol("dt0_theta1")
dt0t0_theta0 = sp.Symbol("dt0t0_theta0")
dt0t1_theta0 = sp.Symbol("dt0t1_theta0")
dt0t0_theta1 = sp.Symbol("dt0t0_theta1")
dt1t0_theta0 = sp.Symbol("dt1t0_theta0")

B = sp.Symbol("B") #B, initial beta 
K = sp.Symbol("K")
beta = sp.Symbol("beta")

epsilon_0 = sp.Symbol("epsilon_0") #epsilon_0
E_r = sp.Function("E")(r_0) #radial dependent epsilon, r0^-2
vrterm = - 2 * epsilon_0 * r_1 * r_0**(-3) * dt0_r0**2 + 6 * epsilon_0 * r_1 * r_0**(-4) * dt0_r0**2 - 2 * epsilon_0 * r_0**(-3) * dt0_r0 * dt0_r1 + epsilon_0 * dt1t0_r0 + epsilon_0 * dt0t1_r0  
omegaterm = - 2 * epsilon_0 * r_0**(-3) * dt0_r0 * dt0_theta1 - 2 * epsilon_0 * r_0**(-3) * dt0_r0 * dt0_theta1 - 2 * epsilon_0 * r_0**(-3) * theta_1 * dt0t0_r0 + 6 * epsilon_0 * r_0**(-4) * dt0_r0**2 * theta_1 + epsilon_0 * dt1t0_theta0 + epsilon_0 * dt0t1_theta0

r_exp = r_0 + epsilon_0 * E_r * r_1 #r perturbed expression
vr_exp = epsilon_0 * E_r * dt0_r1 + epsilon_0 * dt1_r0 + dt0_r0 - 2 * epsilon_0 * r_1 * r_0**(-3) * dt0_r0 #v perturbed expression
vrdot_exp = dt0t0_r0 + epsilon_0 * E_r * dt0t0_r1 - 2 * epsilon_0 * r_0**(-3) * dt0_r0 * dt0_r1 + vrterm #vrdot perturbed expression

theta_exp = theta_0 + epsilon_0 * E_r * theta_1 #theta perturbed expression
omega_exp = dt0_theta0 + epsilon_0 * dt1_theta0 + epsilon_0 * E_r * dt0_theta1 - 2 * epsilon_0 * r_0**(-3) * dt0_r0  #omega perturbed expression
omegadot_exp = dt0t0_theta0 + epsilon_0 * E_r * dt0t0_theta1 + omegaterm #omegadot perturbed expression

m_exp = m_0 + epsilon_0 * E_r * m_1
beta = m_0**(-1 / 3) - 1 / 3 * epsilon_0 * E_r * m_1 * m_0**(-4 / 3) 

#radial equation
rad_eq = sp.Eq((1 - B) * r_exp**2 * (vrdot_exp - r_exp * omega_exp**2) , 
               -(1 - beta * B) - 2 * K * epsilon_0 * beta * B * vr_exp) #radial eq of motion

#up to second order expressions
req_lhs = sp.series(rad_eq.lhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) lhs

req_rhs = sp.series(rad_eq.rhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) rhs

rad_eq = (req_lhs - req_rhs).expand()
rad_eq_zeroth_order = rad_eq.coeff(epsilon_0 , 0) #zeroth order total expression
rad_eq_1 = rad_eq.coeff(epsilon_0 , 1) #1 order expression

#angular equation
ang_eq = sp.Eq(r_exp * (1 - B) * (r_exp * omegadot_exp + 2 * vr_exp * omega_exp) , -K * epsilon_0 * B * beta * omega_exp) #angular eq of motion

#second order sols
angeq_lhs = sp.series(ang_eq.lhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) lhs

angeq_rhs = sp.series(ang_eq.rhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) rhs

angeq = (angeq_lhs - angeq_rhs).expand()
angeq_zeroth_order = angeq.coeff(epsilon_0 , 0) #zeroth order total expression
angeq_1 = angeq.coeff(epsilon_0 , 1) #1 expression

# print(angeq_1)

r0 = sp.Function("r0")(t1)
r1 = sp.Function("r1")(t0)
theta1 = sp.Function("theta1")(t0)
omega0 = sp.Function("omega0")(t1)
omega1 = sp.Function("omega1")(t0)
vr1 = sp.Function("vr1")(t0)
dt1_r0 = sp.Function("dt1_r0")(t1)
dt1_theta0 = sp.Function("dt1_theta0")(t1)
dt0t0_theta1 = sp.Function("dt0t0_theta1")(t0)
dt0_r1 = sp.Function("dt0_r1")(t0)
beta = sp.Function("beta")(t1)
dt1_omega0 = sp.Function("dt1_omega0")(t1)
m0 = sp.Function("m0")(t1)
m1 = sp.Function("m1")(t0)

# dt0_r1 = r0**2 * vr1 - dt1_r0 * r0**2
# dt0_theta1 = r0**2 * omega1 - dt1_theta0 * r0**2

dt0_vr1 = - B * m1 / (3 * (1 - B) * r0**2 * m0**(1 / 3)) + 3 * omega0**2 * r1 + 2 * omega0 * r0 * omega1
dt0_omega1 = - B * K * omega0 / ((1 - B) * m0**(1 / 3) * r0**2) - 2 * omega0 * dt0_r1 / r0**3 - 2 * omega0 * dt1_r0 / r0**2 - dt1_omega0

eqs = [
    # sp.Eq(sp.diff(r1 , t0) , dt0_r1) , 
    # sp.Eq(sp.diff(theta1 , t0) , dt0_theta1) ,
    sp.Eq(sp.diff(vr1 , t0) , dt0_vr1) ,
    sp.Eq(sp.diff(omega1 , t0) , dt0_omega1)
]

sol = sp.dsolve(eqs)

# r1_sol = sol[0]
# theta1_sol = sol[1]
# vr1_sol = sol[0]
# omega1_sol = sol[0]

print(sol[-1])
