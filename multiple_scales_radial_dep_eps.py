import sympy as sp

"""Solves algebraically the equations of motion using perturbed parameters in epsilon_0,
two timescales"""

"""Defining the hatted variables"""
t0 = sp.Symbol("t0") #fast time
t1 = sp.Symbol("t1") #slow time
dt0 = sp.Symbol("dt0") #partial differential, fast time
dt1 = sp.Symbol("dt1") #partial differential, slow time

r_0 = sp.Function("r_0")(t1) #r0
r_1 = sp.Function("r_1")(t0 , t1) #r1

v_r0 = sp.Function("v_r0")(t1) #vr0
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
m_1 = sp.Function("m_1")(t1) #m1

dt0_r0 = sp.Symbol("dt0_r0")
dt0_r1 = sp.Symbol("dt0_r1")
dt1_r0 = sp.Symbol("dt1_r0")
dt1_r1 = sp.Symbol("dt1_r1")

dt0t0_r0 = sp.Symbol("dt0t0_r0")
dt0t0_r1 = sp.Symbol("dt0t0_r1")
dt0t1_r0 = sp.Symbol("dt0t1_r0")
dt1t0_r0 = sp.Symbol("dt1t0_r0")
dt0t1_r1 = sp.Symbol("dt0t1_r1")
dt1t0_r1 = sp.Symbol("dt1t0_r1")
dt1t1_r0 = sp.Symbol("dt1t1_r0")

dt0_theta0 = sp.Symbol("dt0_theta0")
dt1_theta0 = sp.Symbol("dt1_theta0")
dt0_theta1 = sp.Symbol("dt0_theta1")
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
K = sp.Symbol("K") #K, drag-to-mass loss ratio
beta = sp.Symbol("beta") #beta

epsilon_0 = sp.Symbol("epsilon_0") #epsilon_0

r_exp = r_0 + epsilon_0 * r_1 #r perturbed expression

vr_exp = dt0_r0 + epsilon_0 * dt0_r1 + epsilon_0 * dt1_r0 #vr perturbed expression

vrdot_exp = (dt0t0_r0 + epsilon_0 * dt0t0_r1 + epsilon_0 * dt1t0_r0 
             + epsilon_0 * dt0t1_r0) #vrdot perturbed expression

theta_exp = theta_0 + epsilon_0 * theta_1 #theta perturbed expression

omega_exp = dt0_theta0 + epsilon_0 * dt1_theta0 + epsilon_0 * dt0_theta1  #omega perturbed expression

omegadot_exp = (dt0t0_theta0 + epsilon_0 * dt0t0_theta1 + epsilon_0 * dt1t0_theta0 
                + epsilon_0 * dt0t1_theta0) #omegadot perturbed expression

beta = beta_0 + epsilon_0 * beta_1 

#Use the following only after zeroth order
r_exp = r_exp.subs({"dt0_r0" : 0 , "dt0t0_r0" : 0 , "dt0t0_theta0" : 0 , "dt0t1_r0" : 0 , "dt1t0_r0" : 0})
vr_exp = vr_exp.subs({"dt0_r0" : 0 , "dt0t0_r0" : 0 , "dt0t0_theta0" : 0 , "dt0t1_r0" : 0 , "dt1t0_r0" : 0})
vrdot_exp = vrdot_exp.subs({"dt0_r0" : 0 , "dt0t0_r0" : 0 , "dt0t0_theta0" : 0 , "dt0t1_r0" : 0 , "dt1t0_r0" : 0})
theta_exp = theta_exp.subs({"dt0_r0" : 0 , "dt0t0_r0" : 0 , "dt0t0_theta0" : 0 , "dt0t1_r0" : 0 , "dt1t0_r0" : 0})
omega_exp = omega_exp.subs({"dt0_r0" : 0 , "dt0t0_r0" : 0 , "dt0t0_theta0" : 0 , "dt0t1_r0" : 0 , "dt1t0_r0" : 0})
omegadot_exp = omegadot_exp.subs({"dt0_r0" : 0 , "dt0t0_r0" : 0 , "dt0t0_theta0" : 0 , "dt0t1_r0" : 0 , "dt1t0_r0" : 0})

"""Radial equation"""
rad_eq = sp.Eq((1 - B) * r_exp**2 * (vrdot_exp - r_exp * omega_exp**2) , 
               -(1 - beta * B) - 2 * K * epsilon_0 * beta * B * vr_exp) 

"""Keeping terms only up to epsilon^1 order"""
req_lhs = sp.series(rad_eq.lhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) lhs
req_rhs = sp.series(rad_eq.rhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) rhs

"""Separating solutions by order in epsilon"""
rad_eq = (req_lhs - req_rhs).expand()
rad_eq_zeroth_order = rad_eq.coeff(epsilon_0 , 0) #zeroth order in epsilon
rad_eq_1 = rad_eq.coeff(epsilon_0 , 1) #first order in epsilon

"""Angular equation"""
ang_eq = sp.Eq(r_exp * (1 - B) * (r_exp * omegadot_exp + 2 * vr_exp * omega_exp) , 
               -K * epsilon_0 * B * beta * omega_exp) 

"""Keeping terms only up to epsilon^1 order"""
angeq_lhs = sp.series(ang_eq.lhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) lhs
angeq_rhs = sp.series(ang_eq.rhs , epsilon_0 , 0 , 2).removeO() #removing O(epsilon^2) rhs

"""Separating solutions by order in epsilon"""
angeq = (angeq_lhs - angeq_rhs).expand()
angeq_zeroth_order = angeq.coeff(epsilon_0 , 0) #zeroth order in epsilon
angeq_1 = angeq.coeff(epsilon_0 , 1) #first order in epsilon

"""Fetching terms that become secular. Defining variables in two dimensions as functions of t0
as that is of most importance. Solving only dt0_vr1 and dt0_omega1 for cleaner expressions.
In solutions the secular terms are those written as ...*Integral(1,t0)"""
r1 = sp.Function("r1")(t0)
dt1_r0 = sp.Function("dt1_r0")(t1)
dt0_r1 = sp.Function("dt0_r1")(t0)
dt1_r1 = sp.Function("dt1_r1")(t1)
dt0t0_r1 = sp.Function("dt0t0_r1")(t0)

vr1 = sp.Function("vr1")(t0)
dt1_vr1 = sp.Function("dt1_vr1")(t1)

theta1 = sp.Function("theta1")(t0)
dt0_theta1 = sp.Function("dt0_theta1")(t0)
dt1_theta1 = sp.Function("dt1_theta1")(t1)

omega1 = sp.Function("omega1")(t0)

dt1_theta0 = sp.Function("dt1_theta0")(t1)
dt0t0_theta1 = sp.Function("dt0t0_theta1")(t0)

dt1_omega0 = sp.Function("dt1_omega0")(t1)

# dt0_r1 = r_0**2 * vr1 - dt1_r0 * r_0**2
# dt0_theta1 = r_0**2 * omega1 - dt1_theta0 * r_0**2
dt0_vr1 = - B * beta_1 / ((1 - B) * r_0**2) + 3 * omega_0**2 * r1 + 2 * omega_0 * r_0 * omega1
dt0_omega1 = (- 2 * omega_0 * dt0_r1 / r_0 - 2 * omega_0 * dt1_r0 / r_0 - dt1_omega0 
            - B * K * beta_0 * omega_0 / ((1 - B) * r_0**2)) 



eqs = [
    # sp.Eq(sp.diff(r1 , t0) , dt0_r1) , 
    # sp.Eq(sp.diff(theta1 , t0) , dt0_theta1) ,
    sp.Eq(sp.diff(vr1 , t0) , dt0_vr1) ,
    sp.Eq(sp.diff(omega1 , t0) , dt0_omega1)
]

sol = sp.dsolve(eqs)

