import sympy as sp
from sympy import symbols

#defining the hatted variables
t0 , t1 , t2 = symbols("t0") , symbols("t1") , symbols("t2")
dt0 , dt1 , dt2 = symbols("dt0") , symbols("dt1") , symbols("dt2")

r0 = sp.Function("r0")(t0 , t1 , t2) 
r1 = sp.Function("r1")(t0 , t1 , t2) 
r2 = sp.Function("r2")(t0 , t1 , t2) 

vr0 = sp.Function("vr0")(t0 , t1 , t2) 
vr1 = sp.Function("vr1")(t0 , t1 , t2) 
vr2 = sp.Function("vr2")(t0 , t1 , t2) 

theta0 = sp.Function("theta0")(t0 , t1 , t2) 
theta1 = sp.Function("theta1")(t0 , t1 , t2) 
theta2 = sp.Function("theta2")(t0 , t1 , t2) 

omega0 = sp.Function("omega0")(t0 , t1 , t2) 
omega1 = sp.Function("omega1")(t0 , t1 , t2) 
omega2 = sp.Function("omega2")(t0 , t1 , t2) 

beta = symbols("beta")
B = symbols("B")
K = symbols("K")

epsilon = symbols("epsilon") #epsilon

def dt(expr):
    return (sp.diff(expr, t0)
          + epsilon * sp.diff(expr, t1)
          + epsilon**2 * sp.diff(expr, t2))

def dt2(expr):
    return dt(dt(expr))

r = r0 + epsilon * r1 + epsilon**2 * r2
vr = dt(r)
theta = theta0 + epsilon * theta1 + epsilon**2 * theta2
omega = dt(theta)

#radial equation
rad_eq = sp.Eq((vr - r * omega**2) * (1 - B) * r**2 , 
               -(1 - B * beta) * (1 - (2 * K * epsilon * vr))) #radial eq of motion

#up to second order expressions
rad_expr = (rad_eq.lhs - rad_eq.rhs).expand()

rad_eq_zeroth_order = rad_expr.coeff(epsilon , 0) #zeroth order total expression
rad_eq_1 = rad_expr.coeff(epsilon , 1) #1 expression
rad_eq_2 = rad_expr.coeff(epsilon , 2) #2 expression

#angular equation
ang_eq = sp.Eq(r * (1 - B) * (r * dt(omega) + 2 * vr * omega) , (1 - B * beta) * K * epsilon * omega) #angular eq of motion

#second order sols
ang_expr = (ang_eq.lhs - ang_eq.rhs).expand()

angeq_zeroth_order = ang_expr.coeff(epsilon , 0) #zeroth order total expression
angeq_1 = ang_expr.coeff(epsilon , 1) #1 expression
angeq_2 = ang_expr.coeff(epsilon , 2) #2 expression

print(angeq_1)
