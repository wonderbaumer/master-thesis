import sympy as sp

t0 = sp.Symbol("t0") #fast time
t1 = sp.Symbol("t1") #slow time, symbol only

B = sp.Symbol("B") #B, initial beta 
K = sp.Symbol("K")

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

beta = sp.Function("beta")(t1)
dt1_omega0 = sp.Function("dt1_omega0")(t1)
m0 = sp.Function("m0")(t1)
m1 = sp.Function("m1")(t1)
m2 = sp.Function("m2")(t1)

# dt0_r1 = r0**2 * vr1 - dt1_r0 * r0**2
# dt0_theta1 = r0**2 * omega1 - dt1_theta0 * r0**2

dt0_vr1 = - B * m1 / (3 * (1 - B) * r0**2 * m0**(1 / 3)) + 3 * omega0**2 * r1 + 2 * omega0 * r0 * omega1
dt0_omega1 = - B * K * omega0 / ((1 - B) * m0**(1 / 3) * r0**2) - 2 * omega0 * dt0_r1 / r0**3 - 2 * omega0 * dt1_r0 / r0**2 - dt1_omega0

dt0_vr2 = (+ 3 * omega0 * r2 / r0**2 + 3 * omega0**2 * r1**2 / r0**3 #- 2 * B * K * vr1 / ((1 - B) * m0**(1 / 3)) 
           + 6 * omega0 * r1 * dt0_theta1 / r0**2 - 4 * omega0 * dt1_r0 + 6 * omega0 * r1 * dt1_theta0 
           + 2 * omega0 * r0 * dt1_theta1 + dt0_theta1**2 / r0 + 2 * r0 * dt0_theta1 * dt1_theta0 
           - 2 * r1 * dt0t0_r1 / r0**3 + 4 * dt1_r0 * dt1_r1 / r0 + r0**3 * dt1_theta0**2 - r0**2 * dt1_vr1 
           + 2 * r0**2 / (1 - B)) #- B * m2 / (3 * (1 - B) * r0**4 * m0**(1 / 3)) ) 

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

print(sol[-1])