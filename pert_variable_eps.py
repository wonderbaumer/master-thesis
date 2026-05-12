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

    def __init__(self , particle , t , beta):

        self.particle = particle
        self.material = particle.material
        self.size = particle.size

        self.r = particle.r
        self.B = particle.B
        self.V = particle.V
        self.T = particle.T
        self.time = t
        
        self.barr = beta
        self.epsilon0 = particle.eps()
        self.delta = self.V / c
        self.K = self.delta / self.epsilon0
        
        self.beta_prime = np.gradient(self.barr , self.time)
        self.t1 = self.time * self.epsilon0

    def C0(self):
        beta = self.barr
        t1 = self.time * self.epsilon0
        
        cst = -4 * self.K * self.B / (1 - self.B)**3
        

        var = beta * (1 - self.B * beta)**2
        init = 1 / cst
        
        b_int = it.cumulative_simpson(var , x = t1 , initial = 0)
        terms = 1.0 + b_int * cst
        
        invalid = np.where(terms <= 0)
        terms[invalid] = 0.0000001
        
        # C0_tot = terms**(1 / 4)
        C0_tot = np.sqrt(1 - self.barr * self.B)

        return C0_tot

    def C0_prime(self):
        beta = self.barr
        c0 = self.C0()

        C0primetot = -self.B * self.K / (1 - self.B)**3 * beta * (1 - self.B * beta)**2 * c0**(-3)

        return C0primetot
    
    def C3(self):
        
        tot = self.B / (1 - self.B) * (2 * self.K - self.beta_prime[0])
        
        return tot
    
    def D0(self):
        c3 = self.C3()
        tot = -2 * self.epsilon0 * c3

        return tot

    def omega0(self):
        beta = self.barr
        coeff0 = self.C0()

        omega0 = ((1 - self.B) / (1 - beta * self.B))**(-2) * coeff0**(-3)

        return omega0

    def rad(self):
        beta = self.barr
        omega0 = self.omega0()

        coeff0 = self.C0()
        coeff3 = self.C3()
        
        r0 = (1 - self.B) / (1 - beta * self.B) * coeff0**2
        r1 = coeff3 * np.sin(omega0 * self.time)

        firstorder = self.epsilon0 * r0**(-2)
        invalid = np.where(r0 <= 0.1)
        firstorder[invalid] = 0.1
        
        rtot = r0 +  firstorder * r1

        return rtot , r0 , r1

    def omega1(self):
        r0 = self.rad()
        omega0 = self.omega()
        c3 = self.C3()
        

        omega1 = -2 * omega0 / r0**3 * c3 * np.sin(omega0 * self.time)

        return omega1 * self.epsilon0 * r0**(-2)

    def theta(self):
        coeff0 = self.C0()
        coeff0_prime = self.C0_prime()
        c3 = self.C3()
        d0 = self.D0()

        omega0 = self.omega()
        _ , r0 , _ = self.rad()

        theta0 = ((1 - self.B) / (1 - self.B * self.barr))**(-2) * coeff0**(-3) * self.time + d0
        theta11 = 2 * omega0 / r0 * c3 * np.cos(omega0 * self.time)
        theta12 = r0**2 * self.time**2 * self.B * (1 - self.barr * self.B) * self.beta_prime / (1 - self.B)**2 * coeff0**(-3) 
        theta13 = 3 * r0**2 * self.time**2 * (1 - self.barr * self.B) * coeff0_prime / (2 * (1 - self.B) * coeff0**4)

        theta1 = theta11 + theta12 + theta13

        thetatot = theta0 + self.epsilon0 * r0**(-2) * theta1
    
        return thetatot

    def vr(self):
        coeff0 = self.C0()
        coeff0_prime = self.C0_prime()
        c3 = self.C3()
        r0 = self.rad()
        

        _ , omega0 , _ = self.omega()

        vr0 = 0
        vr11 = omega0 / r0**2 * c3 * np.cos(omega0 * self.time)
        vr12 = self.B * (1 - self.B) * coeff0**2 * self.beta_prime / (1 - self.B * self.barr)**2 + 2 * (1 - self.B) * coeff0 * coeff0_prime / (1 - self.B * self.barr)
        vr1 = vr11 + vr12
        
        vrtot = vr0 + self.epsilon0 * r0**(-2) * vr1

        return vrtot  

    def K_rcst(self):
        var = self.barr * (1 - self.B * self.barr)**2
        
        b_int = it.cumulative_simpson(var , x = self.t1 , initial = 0)
        
        K_func = self.beta_prime * (1 - self.B) / (2 * self.barr * (1 - self.B * self.barr))


        return K_func

if __name__== "__main__":
    par = dust_properties("silicate" , "slow" , "large")
    res = np.load("Files/rk45_t7_large_silicate_slowsw_1AU_gradient.npz")
    x , y , _ , _ , m , b , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    
    p = perturbed_functions(par , t , b)
    
    c0 = p.C0()
    
    
    rnum = np.sqrt(x**2+y**2)
    om = p.omega0()
    r , r0 , r1 = p.rad()
    print(p.K , p.K_rcst())
    # print(p.K_rcst() , p.K)
    # plt.plot(t , r0)
    # plt.plot(t , rnum)
    # # plt.ylim(0.0 , 1.01)
    # plt.show()
    
    