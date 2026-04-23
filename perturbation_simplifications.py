import sympy as sp
from sympy import symbols , abc , dsolve , Derivative , integrate
from sympy.abc import a , b , c , g , e , d , f , h , i , j , k , l , m , n , o , p , q , r
from sympy.solvers.ode.systems import dsolve_system
from config import rhat0 , betahat0
from sympy.matrices.expressions import MatMul

t0 = symbols("t0")
t1 = symbols("t1")

b = symbols("b")
c0 = sp.Function("c0")(t1)
B = symbols("B")
K = symbols("K")

beta = 1 / (1 - t1 / 3)
totc0beta = beta * (1 - beta * B)**2
totc0int = integrate(totc0beta , t1)
print(totc0int)

r0_exp = ((1 - B) / (1 - beta * B)) * c0**2
dt1_r0 = sp.diff(r0_exp , t1)

omega0_exp = ((1 - B) / (1 - beta * B))**(-2) * c0**(-3)
dt1_omega0 = sp.diff(omega0_exp , t1)


C0_exp = (-4 * B * K / (1 - B)**3 *(3 * beta**4 * B / 4 - 4 * B**3 * beta**3 + 9 * B**2 * beta**2 - 12 * B * beta + 3 * sp.log(beta)) + 1 + 4 * B * K / (1 - B)**3 * (9 * B**2 - 45 * B / 4 - 4 * B**3))**(1 / 4)

dt1_C0 = sp.diff(C0_exp , t1)


"""
#defining the hatted variables
t = sp.Symbol("t") #time
r_0 = sp.Function("r_0")(t) #r0
r_10 = sp.Function("r_10")(t) #r10
r_01 = sp.Function("r_01")(t) #r01
r_11 = sp.Function("r_11")(t) #r11
r_20 = sp.Function("r_20")(t) #r20 
r_02 = sp.Function("r_02")(t) #r02

v_r0 = sp.Function("v_r0")(t) #v0
v_r10 = sp.Function("v_r10")(t) #vr10
v_r01 = sp.Function("v_r01")(t) #vr01
v_r11 = sp.Function("v_r11")(t) #vr11
v_r20 = sp.Function("v_r20")(t) #vr20
v_r02 = sp.Function("v_r02")(t) #vr02

vrdot_0 = sp.Function("vrdot_0")(t) #vdot0
vrdot_10 = sp.Function("vrdot_10")(t) #vrdot10
vrdot_01 = sp.Function("vrdot_01")(t) #vrdot01
vrdot_11 = sp.Function("vrdot_11")(t) #vrdot11
vrdot_20 = sp.Function("vrdot_20")(t) #vrdot20
vrdot_02 = sp.Function("vrdot_02")(t) #vrdot02

theta_0 = sp.Function("theta_0")(t) #theta0
theta_10 = sp.Function("theta_10")(t) #theta10
theta_01 = sp.Function("theta_01")(t) #theta01
theta_11 = sp.Function("theta_11")(t) #theta11
theta_20 = sp.Function("theta_20")(t) #theta20
theta_02 = sp.Function("theta_02")(t) #theta02

omega_0 = sp.Function("omega_0")(t) #omega0
omega_10 = sp.Function("omega_10")(t) #omega10
omega_01 = sp.Function("omega_01")(t) #omega01
omega_11 = sp.Function("omega_11")(t) #omega11
omega_20 = sp.Function("omega_20")(t) #omega20
omega_02 = sp.Function("omega_02")(t) #omega02


omegadot_0 = sp.Function("omegadot_0")(t) #omegadot0
omegadot_10 = sp.Function("omegadot_10")(t) #omegadot10
omegadot_01 = sp.Function("omegadot_01")(t) #omegadot01
omegadot_11 = sp.Function("omegadot_11")(t) #omegadot11
omegadot_20 = sp.Function("omegadot_20")(t) #omegadot20
omegadot_02 = sp.Function("omegadot_02")(t) #omegadot02

beta_0 = sp.Symbol("beta_0") #beta0
beta_10 = sp.Symbol("beta_10") #beta10
beta_01 = sp.Symbol("beta_01") #beta01
beta_11 = sp.Symbol("beta_11") #beta11
beta_20 = sp.Symbol("beta_20") #beta20
beta_02 = sp.Symbol("beta_02") #beta02


B = sp.Symbol("B") #B, initial beta 
V = sp.Symbol("V") #V, initial velocity
c = sp.Symbol("c") #c, speed of light

epsilon = sp.Symbol("epsilon") #epsilon
delta = sp.Symbol("delta") #delta


r = r_0 + epsilon * r_10 + delta * r_01 + epsilon * delta * r_11 + epsilon**2 * r_20 + delta**2 * r_02 #r perturbed expression
vr = v_r0 + epsilon * v_r10 + delta * v_r01 + epsilon * delta * v_r11 + epsilon**2 * v_r20 + delta**2 * v_r02 #v perturbed expression
vr_dot = vrdot_0 + epsilon * vrdot_10 + delta * vrdot_01 + epsilon * delta * vrdot_11 + epsilon**2 * vrdot_20 + delta**2 * vrdot_02 #vdot perturbed expression

beta = beta_0 + epsilon * beta_10 + delta * beta_01 + epsilon * delta * beta_11 + epsilon**2 * beta_20 + delta**2 * beta_02 #beta petrubed expression

theta = theta_0 + epsilon * theta_10 + delta * theta_01 + epsilon * delta * theta_11 + epsilon**2 * theta_20 + delta**2 * theta_02 #theta perturbed expression
omega = omega_0 + epsilon * omega_10 + delta * omega_01 + epsilon * delta * omega_11 + epsilon**2 * omega_20 + delta**2 * omega_02 #omega perturbed expression
omega_dot = omegadot_0 + epsilon * omegadot_10 + delta * omegadot_01 + epsilon * delta * omegadot_11 + epsilon**2 * omegadot_20 + delta**2 * omegadot_02 #omegadot perturbed expression

#radial equation
rad_eq = sp.Eq((vr_dot - r * omega**2) * (1 - B) * r**2 , 
               -1 + beta * B - 2 * beta * B * vr * delta) #radial eq of motion

#up to second order expressions
req_lhs = sp.series(rad_eq.lhs , epsilon , 0 , 3).removeO() #removing O(epsilon^3) lhs
req_lhs = sp.series(req_lhs , delta , 0 , 3).removeO() #removing O(delta^3) lhs

req_rhs = sp.series(rad_eq.rhs , epsilon , 0 , 3).removeO() #removing O(epsilon^3) rhs
req_rhs = sp.series(req_rhs , delta , 0 , 3).removeO() #removing O(delta^3) rhs


rad_eq = (req_lhs - req_rhs).expand()
rad_eq_zeroth_order = rad_eq.coeff(epsilon , 0).coeff(delta , 0) #zeroth order total expression
rad_eq_10 = rad_eq.coeff(epsilon , 1).coeff(delta , 0) #10 expression
rad_eq_01 = rad_eq.coeff(epsilon , 0).coeff(delta , 1) #01 expression
rad_eq_11 = rad_eq.coeff(epsilon , 1).coeff(delta , 1) #11 expression
rad_eq_20 = rad_eq.coeff(epsilon , 2).coeff(delta , 0) #20 expression
rad_eq_02 = rad_eq.coeff(epsilon , 0).coeff(delta , 2) #02 expression

vrdot0_sol = sp.solve(rad_eq_zeroth_order , vrdot_0) #zeroth order, solving for vdot0
vrdot10_sol = sp.solve(rad_eq_10 , vrdot_10) #first order, solving for vdot10
vrdot01_sol = sp.solve(rad_eq_01 , vrdot_01) #first order, solving for vdot01
vrdot11_sol = sp.solve(rad_eq_11 , vrdot_11) #second order, solving for vdot11
vrdot20_sol = sp.solve(rad_eq_20 , vrdot_20) #second order, solving for vdot20
vrdot02_sol = sp.solve(rad_eq_02 , vrdot_02) #second order, solving for vdot02

#angular equation
ang_eq = sp.Eq(r * (1 - B) * (r * omega_dot + 2 * vr * omega) , - B * beta * omega * delta) #angular eq of motion

#second order sols
angeq_lhs = sp.series(ang_eq.lhs , epsilon , 0 , 3).removeO() #removing O(epsilon^3) lhs
angeq_lhs = sp.series(angeq_lhs , delta , 0 , 3).removeO() #removing O(delta^3) lhs

angeq_rhs = sp.series(ang_eq.rhs , epsilon , 0 , 3).removeO() #removing O(epsilon^3) rhs
angeq_rhs = sp.series(angeq_rhs , delta , 0 , 3).removeO() #removing O(delta^3) rhs

angeq = (angeq_lhs - angeq_rhs).expand()
angeq_zeroth_order = angeq.coeff(epsilon , 0).coeff(delta , 0) #zeroth order total expression
angeq_10 = angeq.coeff(epsilon , 1).coeff(delta , 0) #10 expression
angeq_01 = angeq.coeff(epsilon , 0).coeff(delta , 1) #01 expression
angeq_11 = angeq.coeff(epsilon , 1).coeff(delta , 1) #11 expression
angeq_20 = angeq.coeff(epsilon , 2).coeff(delta , 0) #20 expression
angeq_02 = angeq.coeff(epsilon , 0).coeff(delta , 2) #02 expression

omegadot_0 = sp.solve(angeq_zeroth_order , omegadot_0) #zeroth order, omegadot0
omegadot_10 = sp.solve(angeq_10 , omegadot_10) #first order, omegadot10
omegadot_01 = sp.solve(angeq_01 , omegadot_01) #first order, omegadot01
omegadot_11 = sp.solve(angeq_11 , omegadot_11) #second order, omegadot11
omegadot_20 = sp.solve(angeq_20 , omegadot_20) #second order, omegadot20
omegadot_02 = sp.solve(angeq_02 , omegadot_02) #second order, omegadot02

lhs_angeq_0 = angeq_lhs.expand().coeff(epsilon , 0).coeff(delta , 0)
lhs_angeq_10 = angeq_lhs.expand().coeff(epsilon , 1).coeff(delta , 0)
lhs_angeq_01 = angeq_lhs.expand().coeff(epsilon , 0).coeff(delta , 1)

#print(omegadot_10)


vr1 = sp.Function("vr1")
vr2 = sp.Function("vr2")
theta0 = sp.Function("theta0")
theta1 = sp.Function("theta1")
theta2 = sp.Function("theta2")
r0 = sp.Function("r0")
r1 = sp.Function("r1")
r2 = sp.Function("r2")
omega0 = sp.Function("omega0")
omega1 = sp.Function("omega1")
omega2 = sp.Function("omega2")



dt1_r0 = symbols("dt1_r0")
dt0t1_r1 = symbols("dt0t1_r1")
dt1t1_r0 = symbols("dt1t1_r0")
dt0t0_r1 = symbols("dt0t0_r1")
dt1_r1 = symbols("dt1_r1")
dt2_r0 = symbols("dt2_r0")
dt1t0_r1 = symbols("dt1t0_r1")
dt0t2_r0 = symbols("dt0t2_r0")
dt0t0_r2 = symbols("dt0t0_r2")
dt0_r1 = symbols("dt0_r1")
dt0_r2 = symbols("dt0_r2")
dt1_vr1 = symbols("dt1_vr1")
dt0_vr1 = symbols("dt0_vr1")

dt1_theta0 = symbols("dt1_theta0")
dt0_theta1 = symbols("dt0_theta1")
dt1t0_theta1 = symbols("dt1t0_theta1")
dt0t0_theta1 = symbols("dt0t0_theta1")
dt0t1_theta0 = symbols("dt0t1_theta0")
dt1t1_theta0 = symbols("dt1t1_theta0")
dt1_theta1 = symbols("dt1_theta1")
dt2_theta0 = symbols("dt2_theta0")

dt1_omega0 = symbols("dt1_omega0")
dt2_omega0 = symbols("dt2_omega0")
dt1_omega1 = symbols("dt1_omega1")
dt0_omega1 = symbols("dt0_omega1")
K = symbols("K")
C1 , C2 , C3 , B , beta , A , B , C , D = symbols("C1 C2 C3 B beta A B C D")
R = sp.Function("R")(t0)   # stands for cos(omega0*t0) part
S = sp.Function("S")(t0)   # stands for sin(omega0*t0) part






#a=1 / r0, g=2omega0r0 , d=(1-betaB)K/(1-B) * omega0/r0^2 , f = omega1^2r0 , h = 2r0omega0omega1r1 , i = 2omega0omega2r0, j = omega0^2 , k = 6omega0 , l = r0
#m= 2/r0, n = [1-Bbeta]K/[1-B]r0^2, o = 2omega0/r0^2, p = 2/r0 , q = 2omega0/r0 , r = 3omega0^2/r0

eqs = [Derivative(r1(t0) , t0) - vr1(t0) + dt1_r0 , Derivative(theta1(t0) , t0) - omega1(t0) + dt1_theta0 , 
       Derivative(vr1(t0) , t0) - 3 * j * r1(t0) - g * omega1(t0) , Derivative(omega1(t0) , t0) + q * vr1(t0) - d + dt1_omega0]

sol = dsolve(eqs , [r1(t0) , theta1(t0) , vr1(t0) , omega1(t0)])



#Legg inn enkle løsninger av førsteorden her, få med t0,t1 dep
omega0 = symbols("omega0")
r1_sol = C2 * R + C3 * S + C1 / q 
omega1_sol = -3 * C1  - q * (C2 * R + C3 * S) #
theta1_sol = - m * (C2 * S - C3 * R)
vr1_sol = omega0 * (C3 * R - C2 * S) + (7 * B * beta**(4 / 3) * (beta**(-1) - B)**(2 / 3)) / (9 * (1 - B)**(1 / 3))

der_r2 = Derivative(r2(t0) , t0) - vr2(t0) + dt1_r1 + dt2_r0
dr_theta2 = Derivative(theta2(t0) , t0) - omega2(t0) + dt1_theta1 + dt2_theta0
der_vr2 = Derivative(vr2(t0) , t0) - 2 * n * vr1_sol + dt1_vr1 - l * omega1_sol**2 - k * omega1_sol * r1_sol - g * omega2(t0)  - 3 * r2(t0) * j + m * r1_sol * dt0_vr1 - r * r1_sol**2
der_omega2 = Derivative(omega2(t0) , t0) - n * omega1_sol + dt1_omega1 + dt2_omega0 + m * r1_sol * dt0_omega1 + m * r1_sol * dt1_omega0 + p * vr1_sol * omega1_sol + q * vr2(t0) + o * r1_sol * vr1_sol

eqs_second = [der_r2 , dr_theta2 , der_vr2 , der_omega2]
sol2 = dsolve(eqs_second , [r2(t0) , theta2(t0) , vr2(t0) , omega2(t0)])

"""
