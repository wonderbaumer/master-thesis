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
r_2 = sp.Function("r_2")(t0 , t1)

v_r0 = sp.Function("v_r0")(t1) #v0
v_r1 = sp.Function("v_r1")(t0 , t1) #vr1
v_r2 = sp.Function("v_r2")(t0 , t1)

vrdot_0 = sp.Function("vrdot_0")(t1) #vdot0
vrdot_1 = sp.Function("vrdot_1")(t0 , t1) #vrdot1
vrdot_2 = sp.Function("vrdot_2")(t0 , t1)

theta_0 = sp.Function("theta_0")(t0 , t1) #theta0
theta_1 = sp.Function("theta_1")(t0 , t1) #theta1
theta_2 = sp.Function("theta_2")(t0 , t1)

omega_0 = sp.Function("omega_0")(t1) #omega0
omega_1 = sp.Function("omega_1")(t0 , t1) #omega1
omega_2 = sp.Function("omega_2")(t0 , t1)

omegadot_0 = sp.Function("omegadot_0")(t1) #omegadot0
omegadot_1 = sp.Function("omegadot_1")(t0 , t1) #omegadot1
omegadot_2 = sp.Function("omegadot_2")(t0 , t1)

m_0 = sp.Function("m_0")(t1) #m0
m_1 = sp.Function("m_1")(t0 , t1) #m1
m_2 = sp.Function("m_2")(t0 , t1)

dt0_r0 = sp.Symbol("dt0_r0")
dt0_r1 = sp.Symbol("dt0_r1")
dt1_r0 = sp.Symbol("dt1_r0")
dt0_r2 = sp.Symbol("dt0_r2")
dt1_r1 = sp.Symbol("dt1_r1")

dt0t0_r0 = sp.Symbol("dt0t0_r0")
dt0t0_r1 = sp.Symbol("dt0t0_r1")
dt0t1_r0 = sp.Symbol("dt0t1_r0")
dt1t0_r0 = sp.Symbol("dt1t0_r0")
dt0t0_r2 = sp.Symbol("dt0t0_r2")
dt0t1_r1 = sp.Symbol("dt0t1_r1")
dt1t0_r1 = sp.Symbol("dt1t0_r1")
dt1t1_r0 = sp.Symbol("dt1t1_r0")

dt0_theta0 = sp.Symbol("dt0_theta0")
dt1_theta0 = sp.Symbol("dt1_theta0")
dt0_theta1 = sp.Symbol("dt0_theta1")
dt0_theta2 = sp.Symbol("dt0_theta2")
dt1_theta1 = sp.Symbol("dt1_theta1")

dt0t0_theta0 = sp.Symbol("dt0t0_theta0")
dt0t1_theta0 = sp.Symbol("dt0t1_theta0")
dt0t0_theta1 = sp.Symbol("dt0t0_theta1")
dt1t0_theta0 = sp.Symbol("dt1t0_theta0")
dt0t1_theta1 = sp.Symbol("dt0t1_theta1")
dt1t0_theta1 = sp.Symbol("dt1t0_theta1")
dt1t1_theta0 = sp.Symbol("dt1t1_theta0")

beta_0 = sp.Function("beta_0")(t1)
beta_1 = sp.Function("beta_1")(t1)
beta_2 = sp.Function("beta_2")(t1)

B = sp.Symbol("B") #B, initial beta 
K = sp.Symbol("K")
beta = sp.Symbol("beta")

epsilon_0 = sp.Symbol("epsilon_0") #epsilon_0
E_r = sp.Function("E")(r_0) #radial dependent epsilon, r0^-2

r_exp = r_0 + epsilon_0 * E_r * r_1 + epsilon_0**2 * E_r**2 * r_2 #r perturbed expression

vr_exp = (dt0_r0 + epsilon_0 * E_r * dt0_r1 + epsilon_0 * dt1_r0 - 2 * epsilon_0 * r_1 * r_0**(-3) * dt0_r0
          + epsilon_0**2 * r_0**(-4) * dt0_r2 - 2 * epsilon_0**2 * r_1 * r_0**(-3) * dt1_r0 
          - 4 * epsilon_0**2 * r_0**(-5) * r_2 * dt0_r0 + epsilon_0**2 * E_r * dt1_r1) #v perturbed expression

vrdot_exp = (dt0t0_r0 + epsilon_0 * E_r * dt0t0_r1 - 2 * epsilon_0 * r_0**(-3) * dt0_r0 * dt0_r1 
             - 2 * epsilon_0 * r_1 * r_0**(-3) * dt0t0_r0 + 6 * epsilon_0 * r_1 * r_0**(-4) * dt0_r0**2 
             - 2 * epsilon_0 * r_0**(-3) * dt0_r0 * dt0_r1 + epsilon_0 * dt1t0_r0 + epsilon_0 * dt0t1_r0
             + epsilon_0**2 * r_0**(-4) * dt0t0_r2 - 4 * epsilon_0**2 * r_0**(-5) * dt0_r2 * dt0_r0 
             - 4 * epsilon_0**2 * r_0**(-5) * r_2 * dt0t0_r0 - 4 * epsilon_0**2 * r_0**(-5) * dt0_r0 * dt0_r2 
             + 20 * epsilon_0**2 * r_0**(-6) * r_2 * dt0_r0**2 + epsilon_0**2 * E_r * dt0t1_r1 
             - 2 * epsilon_0**2 * r_0**(-3) * dt1_r1 * dt0_r0 - 2 * epsilon_0**2 * r_1 * r_0**(-3) * dt0t1_r0 
             + 6 * epsilon_0**2 * r_1 * r_0**(-4) * dt0t1_r0 - 2 * epsilon_0**2 * r_0**(-3) * dt1_r1 * dt1_r0   
             - 2 * epsilon_0**2 * r_0**(-3) * dt1_r1 * dt1_r0 + epsilon_0**2 * r_0**(-2) * dt1t0_r1 
             - 2 * epsilon_0**2 - r_0**(-3) * dt0_r1 * dt0_r0 - 2 * epsilon_0**2 * r_1 * r_0**(-3) * dt1t0_r0 
             + 6 * epsilon_0**2 * r_1 * r_0**(-4) * dt0_r0 * dt1_r0 + epsilon_0**2 * dt1t1_r0) #vrdot perturbed expression

theta_exp = theta_0 + epsilon_0 * E_r * theta_1 + epsilon_0**2 * E_r**2 * theta_2 #theta perturbed expression

omega_exp = (dt0_theta0 + epsilon_0 * dt1_theta0 + epsilon_0 * E_r * dt0_theta1 - 2 * epsilon_0 * r_0**(-3) * dt0_r0*
             + epsilon_0**2 * r_0**(-4) * dt0_theta2 - 4 * epsilon_0**2 * r_0**(-5) * theta_2 * dt0_r0 + epsilon_0**2 * E_r * dt1_theta1 
             - 2 * epsilon_0**2 * r_0**(-3) * dt1_r0)  #omega perturbed expression

omegadot_exp = (dt0t0_theta0 + epsilon_0 * E_r * dt0t0_theta1 - 2 * epsilon_0 * r_0**(-3) * dt0_r0 * dt0_theta1 
                - 2 * epsilon_0 * r_0**(-3) * dt0_r0 * dt0_theta1 - 2 * epsilon_0 * r_0**(-3) * theta_1 * dt0t0_r0 
                + 6 * epsilon_0 * r_0**(-4) * dt0_r0**2 * theta_1 + epsilon_0 * dt1t0_theta0 + epsilon_0 * dt0t1_theta0
                + epsilon_0**2 * r_0**(-4) * dt0t0_theta0 - 4 * epsilon_0**2 * r_0**(-5) * dt0_theta0 * dt0_r0 
                - 4 * epsilon_0**2 * r_0**(-5) * theta_2 * dt0t0_r0 - 4 * epsilon_0**2 * r_0**(-5) * dt0_r0 * dt0_theta2 
                + 20 * epsilon_0**2 * r_0**(-6) * theta_2 * dt0_r0**2 + epsilon_0**2 * E_r * dt0t1_theta1 
                - 2 * epsilon_0**2 * r_0**(-3)* dt1_theta1 * dt0_r0 - 2 * epsilon_0**2 * r_0**(-3) * dt0t1_r0 
                + 6 * epsilon_0**2 * r_0**(-4) * dt1_r0 * dt0_r0 + epsilon_0**2 * E_r * dt1t0_theta1 
                - 2 * epsilon_0**2 * r_0**(-3) * dt0_theta1 * dt1_r0 - 2 * epsilon_0**2 * r_0**(-3) * theta_1 * dt1t0_r0 
                - 2 * epsilon_0**2 * r_0**(-3) * dt1_theta1 * dt0_r0 + 6 * epsilon_0**2 * r_0**(-4) * theta_1 * dt0_r0 * dt1_r0 
                + epsilon_0**2 * dt1t1_theta0) #omegadot perturbed expression

#Use the following only after zeroth order
r_exp = r_exp.subs({"dt0_r0" : 0 , "dt0t0_r0" : 0 , "dt0t0_theta0" : 0 , "dt0t1_r0" : 0 , "dt1t0_r0" : 0})
vr_exp = vr_exp.subs({"dt0_r0" : 0 , "dt0t0_r0" : 0 , "dt0t0_theta0" : 0 , "dt0t1_r0" : 0 , "dt1t0_r0" : 0})
vrdot_exp = vrdot_exp.subs({"dt0_r0" : 0 , "dt0t0_r0" : 0 , "dt0t0_theta0" : 0 , "dt0t1_r0" : 0 , "dt1t0_r0" : 0})
theta_exp = theta_exp.subs({"dt0_r0" : 0 , "dt0t0_r0" : 0 , "dt0t0_theta0" : 0 , "dt0t1_r0" : 0 , "dt1t0_r0" : 0})
omega_exp = omega_exp.subs({"dt0_r0" : 0 , "dt0t0_r0" : 0 , "dt0t0_theta0" : 0 , "dt0t1_r0" : 0 , "dt1t0_r0" : 0})
omegadot_exp = omegadot_exp.subs({"dt0_r0" : 0 , "dt0t0_r0" : 0 , "dt0t0_theta0" : 0 , "dt0t1_r0" : 0 , "dt1t0_r0" : 0})

# m_exp = m_0 + epsilon_0 * E_r * m_1 + epsilon_0**2 * E_r**2 * m_2
beta = beta_0 + epsilon_0 * E_r * beta_1 + epsilon_0**2 * E_r**2 * beta_2

#radial equation
rad_eq = sp.Eq((1 - B) * r_exp**2 * (vrdot_exp - r_exp * omega_exp**2) , 
               -(1 - beta * B) - 2 * K * epsilon_0 * beta * B * vr_exp) #radial eq of motion

#up to second order expressions
req_lhs = sp.series(rad_eq.lhs , epsilon_0 , 0 , 3).removeO() #removing O(epsilon^2) lhs

req_rhs = sp.series(rad_eq.rhs , epsilon_0 , 0 , 3).removeO() #removing O(epsilon^2) rhs

rad_eq = (req_lhs - req_rhs).expand()
rad_eq_zeroth_order = rad_eq.coeff(epsilon_0 , 0) #zeroth order total expression
rad_eq_1 = rad_eq.coeff(epsilon_0 , 1) #1 order expression
radeq_2 = rad_eq.coeff(epsilon_0 , 2)

#angular equation
ang_eq = sp.Eq(r_exp * (1 - B) * (r_exp * omegadot_exp + 2 * vr_exp * omega_exp) , -K * epsilon_0 * B * beta * omega_exp) #angular eq of motion

#second order sols
angeq_lhs = sp.series(ang_eq.lhs , epsilon_0 , 0 , 3).removeO() #removing O(epsilon^2) lhs

angeq_rhs = sp.series(ang_eq.rhs , epsilon_0 , 0 , 3).removeO() #removing O(epsilon^2) rhs

angeq = (angeq_lhs - angeq_rhs).expand()
angeq_zeroth_order = angeq.coeff(epsilon_0 , 0) #zeroth order total expression
angeq_1 = angeq.coeff(epsilon_0 , 1) #1 expression
angeq_2 = angeq.coeff(epsilon_0 , 2)

print(radeq_2)

r0 = sp.Function("r0")(t1)
r1 = sp.Function("r1")(t0)
r2 = sp.Function("r2")(t0)
dt1_r0 = sp.Function("dt1_r0")(t1)
dt0_r1 = sp.Function("dt0_r1")(t0)
dt1_r1 = sp.Function("dt1_r1")(t1)
dt0t0_r1 = sp.Function("dt0t0_r1")(t0)

dt1_vr1 = sp.Function("dt1_vr1")(t1)

theta1 = sp.Function("theta1")(t0)
dt0_theta1 = sp.Function("dt0_theta1")(t0)
dt1_theta1 = sp.Function("dt1_theta1")(t1)

omega0 = sp.Function("omega0")(t1)
omega1 = sp.Function("omega1")(t0)
vr1 = sp.Function("vr1")(t0)
vr2 = sp.Function("vr2")(t0)

dt1_theta0 = sp.Function("dt1_theta0")(t1)
dt0t0_theta1 = sp.Function("dt0t0_theta1")(t0)

beta0 = sp.Function("beta0")(t1)
beta1 = sp.Function("beta1")(t1)
beta2 = sp.Function("beta2")(t1)

dt1_omega0 = sp.Function("dt1_omega0")(t1)
m0 = sp.Function("m0")(t1)
m1 = sp.Function("m1")(t1)
m2 = sp.Function("m2")(t1)

# dt0_r1 = r0**2 * vr1 - dt1_r0 * r0**2
# dt0_theta1 = r0**2 * omega1 - dt1_theta0 * r0**2

# dt0_vr1 = - B * m1 / (3 * (1 - B) * r0**2 * m0**(1 / 3)) + 3 * omega0**2 * r1 + 2 * omega0 * r0 * omega1
# dt0_omega1 = - 2 * omega0 * dt0_r1 / r0**3 #- 2 * omega0 * dt1_r0 / r0**2 - dt1_omega0 - B * K * omega0 / ((1 - B) * m0**(1 / 3) * r0**2) 

dt0_vr1 = sp.Function("dt0_vr1")(t0)
dt0_omega1 = sp.Function("dt0_omega1")(t0)

dt0_vr2 = (+ 3 * omega0 * r2 / r0**2 + 3 * omega0**2 * r1**2 / r0**3 - 2 * B * K * vr1 * beta0 / (1 - B) 
           + 6 * omega0 * r1 * dt0_theta1 / r0**2 - 4 * omega0 * dt1_r0 + 6 * omega0 * r1 * dt1_theta0 
           + 2 * omega0 * r0 * dt1_theta1 + dt0_theta1**2 / r0 + 2 * r0 * dt0_theta1 * dt1_theta0 
           - 2 * r1 * dt0t0_r1 / r0**3 + 4 * dt1_r0 * dt1_r1 / r0 + r0**3 * dt1_theta0**2 - r0**2 * dt1_vr1
             - B * beta2 / ((1 - B) * r0**4)+ 2 * r0**2 / (1 - B))

eqs = [
    # sp.Eq(sp.diff(r1 , t0) , dt0_r1) , 
    # sp.Eq(sp.diff(theta1 , t0) , dt0_theta1) ,
    sp.Eq(sp.diff(vr2 , t0) , dt0_vr2) ,
    # sp.Eq(sp.diff(omega1 , t0) , dt0_omega1)
]

sol = sp.dsolve(eqs)

# r1_sol = sol[0]
# theta1_sol = sol[1]
# vr1_sol = sol[0]
# omega1_sol = sol[0]

# print(sol[-1])
