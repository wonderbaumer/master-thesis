import numpy as np
from config import sil_betaval , car_betaval_bound , init_vals , material_files , m_s , R , size_to_mass , car_size_bound , sil_size
from scipy.constants import c , G
from forces_scaled import inter_func
from dust_properties import dust_properties
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator as pchip
from numpy.lib.stride_tricks import sliding_window_view
from config import t6 , t5
import scipy.integrate as it
import matplotlib.pyplot as plt

class perturbed_functions():

    def __init__(self , particle , t , beta , find_k = False):

        self.particle = particle
        self.material = particle.material
        self.size = particle.size

        self.r = particle.r
        self.B = particle.B
        self.V = particle.V
        self.T = particle.T
        self.time = t
        
        self.barr = beta
        self.epsilon = particle.eps()
        self.delta = self.V / c
        self.K = self.delta / self.epsilon
        self.K_theoretical = self.K_calc()

    
    def K_calc(self):
        
        tot = self.barr * (1 - self.B)**3 / (6 * (1 - self.barr * self.B))

        return tot
    
    def betahat_analytical(self):
        return 1 / (1 - self.time * self.epsilon / 3)
    
    ###WITH DRAG###
    def C0(self , k):
        # beta = self.betahat_analytical()
        beta = self.barr
        t1 = self.time * self.epsilon
        
        cst = -4 * k * self.B / (1 - self.B)**3
        

        var = beta * (1 - self.B * beta)**2
        init = 1 / cst
        
        b_int = it.cumulative_simpson(var , x = t1 , initial = 0)
        terms = 1.0 + b_int * cst

        """
        z = np.log(-1 / beta + 0j)
        first = -3 * z.real 
        second = -(-27 * self.B**2 - 36 * self.B * t1 + 108 * self.B) / (2 * t1**2 - 12 * t1 + 18)

        c = 1 - 2 * self.B**2 * k / (9 * (1 - self.B)**3) * (108 - 27 * self.B)
        terms_man= cst * (first + second) + c
        """

        invalid = np.where(terms <= 0)
        terms[invalid] = 0.0000001
        
        C0_tot = terms**(1 / 4)

        return C0_tot

    def C0_prime(self , k):
        # beta = self.betahat_analytical()
        beta = self.barr
        c0 = self.C0(k)

        C0primetot = -self.B * k / (1 - self.B)**3 * beta * (1 - self.B * beta)**2 * c0**(-3)

        return C0primetot
    
    def C3(self , k):
        tot = self.B * (24 * k - 1) / (3 * (1 - self.B))

        return tot
    
    def D0(self , c3):
        tot = -2 * self.epsilon * c3

        return tot

    def omega(self , k):
        beta = self.barr
        coeff0 = self.C0(k)
        c3 = self.C3(k)

        omega0 = ((1 - self.B) / (1 - beta * self.B))**(-2) * coeff0**(-3)
        omega1 = -2 * omega0 / (((1 - self.B) / (1 - beta * self.B)) * coeff0**2) * c3 * np.sin(omega0 * self.time)

        omegatot = omega0 + self.epsilon * omega1

        return omegatot , omega0 , omega1

    def rad(self , k):
        # beta = self.betahat_analytical()
        beta = self.barr
        _ , omega0 , _ = self.omega(k)

        coeff0 = self.C0(k)
        coeff3 = self.C3(k)
        
        r0 = (1 - self.B) / (1 - beta * self.B) * coeff0**2
        r1 = coeff3 * np.sin(omega0 * self.time)
        
        rtot = r0 + self.epsilon * r1

        return rtot , r0 , r1

    def theta(self , k):
        coeff0 = self.C0(k)
        coeff0_prime = self.C0_prime(k)
        c3 = self.C3(k)
        d0 = self.D0(c3)

        _ , omega0 , _ = self.omega(k)
        _ , r0 , _ = self.rad(k)

        theta0 = ((1 - self.B) / (1 - self.B * self.barr))**(-2) * coeff0**(-3) * self.time + d0
        theta11 = 2 * omega0 / r0 * c3 * np.cos(omega0 * self.time)
        theta12 = (self.time**2 / 2 * coeff0**(-3) * (1 - self.barr * self.B) / (1 - self.B)**2) * (2 * self.B * self.barr**2 / 3 + 3 * (1 - self.barr * self.B) * coeff0**(-1) * coeff0_prime)

        theta1 = theta11 + theta12

        thetatot = theta0 + self.epsilon * theta1
    
        return thetatot

    def vr(self , k):
        coeff0 = self.C0(k)
        coeff0_prime = self.C0_prime(k)
        c3 = self.C3(k)
        

        _ , omega0 , _ = self.omega()

        vr0 = 0
        vr11 = omega0 * c3 * np.cos(omega0 * self.time)
        vr12 = ((1 - self.B) / (1 - self.B * self.barr)) * coeff0 * (self.barr**2 * self.B * coeff0 / (3 * (1 - self.B)) + 2 * coeff0_prime)
        vr1 = vr11 + vr12

        vrtot = vr0 + self.epsilon * vr1

        return vrtot    

if __name__== "__main__":
    par = dust_properties("carbon" , "slow" , "large")
    res = np.load("Files/rk45_t6_large_carbon_slowsw.npz")
    x , y , _ , _ , m , b , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    
    rnum = np.sqrt(x**2+y**2)
    p = perturbed_functions(par , t , b , find_k = False)
    print(p.K_theoretical , par.K)
    # c0 = p.C0(p.K)
    
    
    
    # om , _ , _ = p.omega(p.K)
    # r , _ , _ = p.rad(p.K)


    
    

    
    
    
    
    
   