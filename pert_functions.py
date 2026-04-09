import numpy as np
from config import sil_beta , car_beta , init_vals , material_files
from scipy.constants import c
from forces_scaled import inter_func
from dust_properties import dust_properties
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator as pchip

class perturbed_functions():

    def __init__(self , particle , t , beta):

        self.particle = particle
        self.material = particle.material
        self.size = particle.size

        self.r = particle.r
        self.B = particle.B
        self.V = particle.V
        
        self.time = t
        self.barr = beta
        self.epsilon = particle.eps()
        self.K = self.V / c / self.epsilon
        self.coeff3 = self.C3()
        self.d0 = self.D0()
        #self.beta_func = pchip(self.time , self.barr , extrapolate = False)
        # self.C0_0 = 1.0
        # self.c0 = None
        # self.c0_func = None
        # self.beta_vals = None
        #print(self.epsilon / 2.6031255888589655) #1.9279315785095778e-06)
        
    
    def betahat_analytical(self):
        bvals = self.epsilon * self.time / 2.6031255888589655 + 0.9996231879878003
        return bvals
    
    """
    #calculates up to first order of beta hat from analytical solution
    def betahat_analytical(self):
        ""input: t (array), scaled time to calculate for

        returns: bvals (array), contains all orders as well as their sum""
    
        bvals = 1 / (1 -  self.epsilon * self.time / 3) #total expression
        
        return bvals
    
    
    def du(self , t , u):
        #c0 = c0[0]
        #c0_safe = max(c0, 1e-12)
        beta = self.beta_func(t)
        tot = np.array([-4 * self.B * self.K * beta * (1 - beta * self.B)**4 / ((1 - self.B)**3)])

        return tot
    
    def solve_c0(self):

        ttot = (self.time[0] , self.time[-1])

        sol = solve_ivp(self.du , ttot , [self.C0_0**4] , method = "Radau" , t_eval = self.time , rtol=1e-9 , atol=1e-12)

        self.u = sol.y[0]
        self.beta_vals = self.beta_func(self.time)
        self.c0 = self.u**(1 / 4)

        return self.c0
    
    
    ###WITH DRAG###
    def C0(self , beta_func):
        #beta_func = self.betahat_analytical()
        first = -4 * self.B * self.K / (1 - self.B)**3
        second = (3 * beta_func**4 * self.B / 4 - 4 * self.B**3 * beta_func**3 + 9 * self.B**2 * beta_func**2 - 12 * self.B * beta_func + 3 * np.log(beta_func))
        third = 1 + (4 * self.B * self.K / (1 - self.B)**3) * (9 * self.B**2 - 45 * self.B / 4 - 4 * self.B**3)

        terms = first * second + third
        terms = np.where(terms <= 0 , 0.0000001 , terms)

        tot = terms**(1 / 4)

        return tot
    """
    def C0(self , b_func):
        beta0 = 0.9996231879878003
        a = 2.6031255888589655
        first = -4 * self.B * self.K / (1 - self.B)**3
        second = (a * self.B * b_func**6 / 6 - 4 * self.B**3 * a * b_func**5 / 5 + 6 * a * self.B**2 * b_func**4 / 4 - 4 * a * self.B * b_func**3 / 3 + a * b_func**2 / 2)
        third = 1 + (4 * self.B * self.K / (1 - self.B)**3) * (a * self.B * beta0**6 / 6  - 4 * a * self.B**3 * beta0**5 / 5 + 6 * a / 4 * self.B**2 * beta0**4 - 4 * a / 3 * self.B * beta0**3 + a / 2 * beta0**2)

        terms = first * second + third
        #terms = np.where(terms <= 0 , 0.0000001 , terms)

        tot = terms**(1 / 4)

        return tot
    """
    def C0_prime(self , b_func):

        numerator = -self.B * self.K * (-4 * self.B**3 * b_func**4 + 6 * self.B**2 * b_func**3 - 4 * self.B * b_func**2 + self.B * b_func**5 + b_func)
        denom1 = 4 * self.B * self.K * (-4 * self.B**3 + 9 * self.B**2 - 45 * self.B / 4)
        denom2 = -4 * self.B * self.K * (-4 * self.B**3 * b_func**3 + 9 * self.B**2 * b_func**2 - 12 * self.B * b_func + 3 * self.B * b_func**4 / 4 + 3 * np.log(b_func))
        denom3 = (1 - self.B)**(-3)
        tot_denom = (denom1 + denom2 + denom3)**(3 / 4)

        tot = numerator / tot_denom
        
        return tot
    """
    def C3(self):
        first = -self.B / (3 * (1 - self.B))
        second = -2 * self.B * self.K * (-4 * self.B**3 + 6 * self.B**2 - 3 * self.B + 1) / (1 - self.B)**3

        tot = first + second

        return tot
    
    def D0(self):
        tot = -2 * self.epsilon * self.coeff3

        return tot

    def omega(self):
        b_func = self.betahat_analytical()
        coeff0 = self.C0(b_func)

        omega0 = ((1 - self.B) / (1 - b_func * self.B))**(-2) * coeff0**(-3)
        omega1 = -2 * omega0 / (((1 - self.B) / (1 - b_func * self.B)) * coeff0**2) * self.coeff3 * np.sin(omega0 * self.time)

        omegatot = omega0 + self.epsilon * omega1

        return omegatot , omega0 , omega1

    def rad(self):
        _ , omega0 , _ = self.omega()

        b_func = self.betahat_analytical()
        coeff0 = self.C0(b_func)

        r0 = (1 - self.B) / (1 - b_func * self.B) * coeff0**2
        r1 = self.coeff3 * np.sin(omega0 * self.time)

        rtot = r0 + self.epsilon * r1

        return rtot , r0 , r1

    def theta(self):
        b_func = self.betahat_analytical()
        coeff0 = self.C0(b_func)
        coeff0_prime = 1.0 #self.C0_prime(b_func)

        _ , omega0 , _ = self.omega()
        _ , r0 , _ = self.rad()

        theta0 = ((1 - self.B) / (1 - self.B * b_func))**(-2) * coeff0**(-3) * self.time + self.d0
        theta11 = 2 / r0 * self.coeff3 * np.cos(omega0 * self.time)
        theta12 = -self.time**2 / 2 * coeff0**(-3) * ((1 - b_func * self.B) / (1 - self.B)**2) * (-2 * self.B * b_func**2 / 3 - 3 * (1 - b_func * self.B) * coeff0**(-1) * coeff0_prime)

        theta1 = theta11 + theta12

        thetatot = theta0 + self.epsilon * theta1
    
        return thetatot

    def vr(self):
        b_func = self.betahat_analytical()
        coeff0 = self.C0(b_func)
        #coeff0_prime = self.C0_prime(b_func)
        coeff0_prime = 1.0

        _ , omega0 , _ = self.omega()

        vr0 = 0
        vr11 = omega0 * self.coeff3 * np.cos(omega0 * self.time)
        vr12 = ((1 - self.B) / (1 - self.B * b_func)) * coeff0 * (b_func**2 * self.B * coeff0 / (3 * (1 - self.B)) + 2 * coeff0_prime)
        vr1 = vr11 + vr12

        vrtot = vr0 + self.epsilon * vr1

        return vrtot


if __name__== "__main__":
    par = dust_properties("carbon" , "slow" , "all" , "large")
    res = np.load("Files/rk45_t6_large_carbon_slowsw.npz")
    _ , _ , _ , _ , _ , b , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]

    p = perturbed_functions(par , t , b)
    
    beta = p.betahat_analytical()
    c0 = p.C0(beta)
    print(f"r:" ,p.rad() , "beta:" ,beta ,"c0:", c0)
    
    omegatot , rtot , theta , vr , c0 , beta = p.omega() , p.rad() , p.theta() , p.vr() , c0 , beta
    omega , _ , _ = omegatot
    r , _ , _ = rtot
    np.savez("Files/pert_t6_large_carbon_slowsw.npz" , omega = omega , r = r , theta = theta , vr = vr , c0 = c0 , b = beta , t = t)

    