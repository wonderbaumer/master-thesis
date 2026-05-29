import sympy as sp

t0, t1 = sp.symbols('t0 t1')
B, K = sp.symbols('B K')

# --- slow variables (PARAMETERS in t0-system) ---
r0 = sp.Function('r0')(t1)
omega0 = sp.Function('omega0')(t1)
m0 = sp.Function('m0')(t1)

# derivatives in slow time
dr0_t1 = sp.diff(r0, t1)
domega0_t1 = sp.diff(omega0, t1)

# --- fast-time unknowns ONLY depend on t0 ---
r1 = sp.Function('r1')(t0)
theta1 = sp.Function('theta1')(t0)
vr1 = sp.Function('vr1')(t0)
omega1 = sp.Function('omega1')(t0)

# --- useful structures ---
theta_term = sp.diff(theta1**2, t0, 2)

# --- RHS examples ---
dt0_r1 = r0**2 * vr1 - r0**2 * dr0_t1

dt0_theta1 = r0**2 * omega1  # (plus slow forcing if needed)

dt0_vr1 = (
    -B * m0**(-sp.Rational(1,3)) / (3*(1-B)*r0**2)
    + 3*omega0**2*r1
    + 2*omega0*r0*omega1
)

dt0_omega1 = (
    -B*K*omega0 / ((1-B)*m0**(sp.Rational(1,3))*r0**2)
    - 2*vr1*omega0 / r0**2
    + theta_term / r0**3
    - domega0_t1
)

# --- system in t0 only ---
eqs = [
    sp.Eq(sp.diff(r1, t0), dt0_r1),
    sp.Eq(sp.diff(theta1, t0), dt0_theta1),
    sp.Eq(sp.diff(vr1, t0), dt0_vr1),
    sp.Eq(sp.diff(omega1, t0), dt0_omega1),
]

sol = sp.dsolve(eqs)
print(sol[-1])