import sympy as sp

"""File only used for algebraically differentiating and integrating expressions"""

t0 = sp.Symbol("t0")
t1 = sp.Symbol("t1")

B = sp.Symbol("B")
K = sp.Symbol("K")

c0 = sp.Function("c0")(t1)
d0 = sp.Function("d0")(t1)
m0 = sp.Function("m0")(t1)
m1 = sp.Function("m1")(t1)
r0 = sp.Function("r0")(t1)
beta0 = sp.Function("beta0")(t1)

theta0 = ((1 - B) / (1 - beta0 * B))**(-2) * c0**(-3) * t0 + d0
dt1_theta0 = sp.diff(theta0 , t1)

r0_exp = ((1 - B) / (1 - beta0 * B)) * c0**2
dt1_r0 = sp.diff(r0_exp , t1)

w0_exp = ((1 - B) / (1 - beta0 * B))**(-2) * c0**(-3)
dt1_w0 = sp.diff(w0_exp , t1)

r1_expr = B * m1 / (3 * (1 - B) * r0**2 * m0**(1/3))
dt1_r1exp = sp.diff(r1_expr , t1)

c0 = integrate(m0**(-1 / 3) , t1)











