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
    
    # def C3(self):

    #     tot = self.B / (1 - self.B) * (2 * self.K - self.beta_prime[0])

        
    #     return tot

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
        
        C0_tot = terms**(1 / 4)

        return C0_tot

    def C0_prime(self):
        beta = self.barr
        c0 = self.C0()

        C0primetot = -self.B * self.K / (1 - self.B)**3 * beta * (1 - self.B * beta)**2 * c0**(-3)

        return C0primetot
    

    def omega0(self):
        beta = self.barr
        coeff0 = self.C0()

        omega0 = ((1 - self.B) / (1 - beta * self.B))**(-2) * coeff0**(-3)

        return omega0

    def rad(self):
        beta = self.barr
        omega0 = self.omega0()

        coeff0 = self.C0()
        # coeff3 = self.C3()
        
        r0 = (1 - self.B) / (1 - beta * self.B) * coeff0**2
        # r1 = coeff3 * np.sin(omega0 * self.time)

        firstorder = self.epsilon0 * r0**(-2)
        invalid = np.where(r0 <= 0.1)
        firstorder[invalid] = 0.1
        
        rtot = r0 #+ firstorder * r1

        return rtot 

    def theta(self):
        coeff0 = self.C0()
        coeff0_prime = self.C0_prime()
        d0 = 0.0 #placeholder

        omega0 = self.omega()
        _ , r0 , _ = self.rad()

        theta0 = ((1 - self.B) / (1 - self.B * self.barr))**(-2) * coeff0**(-3) * self.time + d0
        # theta11 = 2 * omega0 / r0 * c3 * np.cos(omega0 * self.time)
        # theta12 = r0**2 * self.time**2 * self.B * (1 - self.barr * self.B) * self.beta_prime / (1 - self.B)**2 * coeff0**(-3) 
        # theta13 = 3 * r0**2 * self.time**2 * (1 - self.barr * self.B) * coeff0_prime / (2 * (1 - self.B) * coeff0**4)

        # theta1 = theta11 + theta12 + theta13

        thetatot = theta0 #+ self.epsilon0 * r0**(-2) * theta1
    
        return thetatot

    def vr(self):
        coeff0 = self.C0()
        coeff0_prime = self.C0_prime()
        r0 = self.rad()
        

        omega0 = self.omega0()

        vr0 = 0
        # vr11 = omega0 / r0**2 * c3 * np.cos(omega0 * self.time)
        vr12 = self.B * (1 - self.B) * coeff0**2 * self.beta_prime / (1 - self.B * self.barr)**2 + 2 * (1 - self.B) * coeff0 * coeff0_prime / (1 - self.B * self.barr)
        # vr1 = vr11 + vr12
        
        vrtot = vr0 + self.epsilon0 * r0**(-2) * vr12

        return vrtot  

    def K_rcst(self):
        var = self.barr * (1 - self.B * self.barr)**2
        
        b_int = it.cumulative_simpson(var , x = self.t1 , initial = 0)
        
        K_func = self.beta_prime * (1 - self.B) / (2 * self.barr * (1 - self.B * self.barr))


        return K_func


def solve_system(t, B, K , epsilon):

    def rhs(t, y):
        beta0 = y[0]
        C0 = y[1]

        m0 = 1.0
        beta0 = m0**(-1 / 3)

        r0 = C0**2 * (1 - B) / (1 - B * beta0)
        dbeta0_dt = beta0**2 / (3 * r0**2)

        dC0_dt = -B * K * beta0 * (1 - beta0 * B)**2 / (1 - B)**3 * C0**(-3)

        return [dbeta0_dt, dC0_dt]

    y0 = [1.0 , 1.0]

    sol = solve_ivp(rhs, (t[0], t[-1]), y0, t_eval=t, method="RK45", rtol=1e-9 , atol=1e-12 )

    beta0 = sol.y[0]
    C0 = sol.y[1]


    r0 = C0**2 * (1 - B) / (1 - B * beta0)
 
    omega0 = C0**(-3) * ((1 - B) / (1 - B * beta0))**(-2)
    

    return r0 , omega0 , beta0 , sol.t

if __name__== "__main__":
    par = dust_properties("silicate" , "slow" , "large")
    res = np.load("Files/rk45_t6_large_silicate_slowsw.npz")
    x , y , _ , _ , m , b , t  = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    
    p = perturbed_functions(par , t , b)
    
    c0 = p.C0()
    
    
    rnum = np.sqrt(x**2+y**2)
    om = p.omega0()
    r = p.rad()
    
    rad , omega0 , beta0 , time = solve_system(p.epsilon0 * t , p.B , p.K , p.epsilon0)
   
    mask = rad >= 0.1
    
    trad = t[:len(rad)]
    xrad = x
    yrad = y
    mrad = m
    rrad = np.sqrt(x**2+y**2)[:len(rad)]
    brad = b[:len(rad)]
    
    plt.plot(trad , rrad , label = "num")
    plt.plot(trad  , rad  , label = "pert" , linestyle = "--")
    plt.legend()
    plt.show()
    
    
    