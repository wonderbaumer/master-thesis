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
        self.K = self.V / (c * self.epsilon)
        # self.coeff3 = self.C3()
        # self.d0 = self.D0()
        self.find_k = find_k
        if self.find_k:
            self.kb_comb()
            #self.K_variation()
    
    def K_calc(self):
        
        tot = (1 - self.B)**5 / (6 * self.barr * (1 - self.barr * self.B)**7)

        return tot
    
    ###WITH DRAG###
    def C0(self , k):
        cst = -4 * self.B * k / (1 - self.B)**3
        
        var = self.barr * (1 - self.B * self.barr)**4
        b_int = it.cumulative_trapezoid(var , self.time , initial = 0)
        terms = 1.0 + cst * b_int * self.epsilon
        invalid = np.where(terms <= 0)
        terms[invalid] = 0.0000001

        C0_tot = terms**(1 / 4)

        return C0_tot
    
    def kb_comb(self):
        b_func_test = self.betahat_analytical(self.epsilon)[::10]
        
        if self.material == "silicate":
           Bvals = sil_betaval
           r_vals = sil_size * 1e-6
           m_vals = size_to_mass(r_vals , "silicate")
        
        if self.material == "carbon":
            Bvals = car_betaval_bound
            r_vals = car_size_bound * 1e-6
            m_vals = size_to_mass(r_vals , "carbon")

        B0 = self.particle.B
        T0 = self.particle.T
        r0 = self.particle.r
        m0 = self.particle.m0

        b_arr = np.zeros(len(Bvals))
        k_vals = np.linspace(0 , 100 , len(b_func_test))
        k_arr = np.zeros((len(Bvals) , len(b_func_test)))
        r0_arr = np.zeros((len(Bvals) , len(b_func_test)))
        #betaval = np.zeros(len(b_func_test))
        # T_arr = np.zeros(len(Bvals))
        # V_arr = np.zeros(len(Bvals))
        # eps_arr = np.zeros(len(Bvals))
        
        for idx , (i , r_i, m_i) in enumerate(zip(Bvals , r_vals, m_vals)):
            T = round(np.sqrt(R**3 / (G * m_s * (1 - i))))

            self.particle.B = i
            self.particle.r = r_i
            self.particle.m0 = m_i
            self.particle.T = T
            epsilon = self.particle.eps()
            b_func = self.betahat_analytical(epsilon)[::10]

            V = np.sqrt((G * m_s * (1 - i)) / R)
            #K = V / (c * epsilon)
            K = np.linspace(0 , 100 , len(Bvals))
            coeff0 = self.C0(b_func , i , k_vals)
        
            r0 = (1 - i) / (1 - b_func * i) * coeff0**2
        
            b_arr[idx] = i
            k_arr[idx , :] = k_vals
            r0_arr[idx , :] = r0
            #betaval[idx] = b_func
            # T_arr[idx] = i
            # V_arr[idx] = i
            # eps_arr[idx] = i

        return b_arr , k_arr , r0_arr #, betaval #, T_arr , V_arr , eps_arr

    def K_variation(self):
        b_func_test = self.betahat_analytical(self.epsilon)[::10]
        
        if self.material == "silicate":
           Bvals = sil_betaval
           r_vals = sil_size * 1e-6
           m_vals = size_to_mass(r_vals , "silicate")
        
        if self.material == "carbon":
            Bvals = car_betaval_bound
            r_vals = car_size_bound * 1e-6
            m_vals = size_to_mass(r_vals , "carbon")

        B0 = self.particle.B
        T0 = self.particle.T
        r0 = self.particle.r
        m0 = self.particle.m0

        b_arr = np.zeros(len(Bvals))
        k_arr = np.zeros((len(Bvals) , len(b_func_test)))
        betaval = np.zeros(len(b_func_test))
        eps_arr = np.zeros(len(Bvals))

        for idx , (i, r_i, m_i) in enumerate(zip(Bvals, r_vals, m_vals)):
            T = round(np.sqrt(R**3 / (G * m_s * (1 - i))))

            self.particle.B = i
            self.particle.r = r_i
            self.particle.m0 = m_i
            self.particle.T = T
            epsilon = self.particle.eps()
            b_func = self.betahat_analytical(epsilon)[::10]
            b_func = np.where(b_func <= 0 , 0.0000001 , b_func)
            num = (((1 - i * b_func) / (1 - i))**2 - 1) * (1 - i)**3
            denom = -3 * b_func**4 * i / 4 + 4 * i**3 * b_func**3 - 9 * i**2 * b_func**2 + 12 * i * b_func - 3 * np.log(b_func) - 4 * i**3 + 9 * i**2 - 45 * i / 4

            K = num / (4 * i * denom)
                

            b_arr[idx] = i
            k_arr[idx , :] = K
            betaval[:] = b_func
            eps_arr[idx] = epsilon

        return b_arr , k_arr , betaval , eps_arr

    def C0_prime(self , k):

        cst = -4 * self.B * k / (1 - self.B)**3

        
        var = self.barr * (1 - self.B * self.barr)**4
        b_int = it.cumulative_trapezoid(var , self.time , initial = 0)
        terms = 1.0 + cst * b_int * self.epsilon
        invalid = np.where(terms <= 0)
        terms[invalid] = 0.0000001

        first = terms**(-3 / 4)

        second = cst * self.barr * (1 - self.B * self.barr)**4
        tot = 1 / 4 * first * second

        return tot
    
    def C3(self , k):
        first = -self.B / (3 * (1 - self.B))
        second = -2 * self.B * k * (-4 * self.B**3 + 6 * self.B**2 - 3 * self.B + 1) / (1 - self.B)**3

        tot = first + second

        return tot
    
    def D0(self , c3):
        tot = -2 * self.epsilon * c3

        return tot

    def omega(self):
        coeff0 = self.C0(self.K)
        c3 = self.C3(self.K)

        omega0 = ((1 - self.B) / (1 - self.barr * self.B))**(-2) * coeff0**(-3)
        omega1 = -2 * omega0 / (((1 - self.B) / (1 - self.barr * self.B)) * coeff0**2) * c3 * np.sin(omega0 * self.time)

        omegatot = omega0 + self.epsilon * omega1

        return omegatot , omega0 , omega1

    def rad(self):
        _ , omega0 , _ = self.omega()

        coeff0 = self.C0(self.K)
        coeff3 = self.C3(self.K)

        r0 = (1 - self.B) / (1 - self.barr * self.B) * coeff0**2
        r1 = coeff3 * np.sin(omega0 * self.time)
        
        rtot = r0 + self.epsilon * r1

        return rtot , r0 , r1

    def theta(self):
        coeff0 = self.C0(self.K)
        coeff0_prime = self.C0_prime(self.K)
        c3 = self.C3(self.K)
        d0 = self.D0(c3)

        _ , omega0 , _ = self.omega()
        _ , r0 , _ = self.rad()

        theta0 = ((1 - self.B) / (1 - self.B * self.barr))**(-2) * coeff0**(-3) * self.time + d0
        theta11 = 2 / r0 * c3 * np.cos(omega0 * self.time)
        theta12 = -self.time**2 / 2 * coeff0**(-3) * ((1 - self.barr * self.B) / (1 - self.B)**2) * (-2 * self.B * self.barr**2 / 3 - 3 * (1 - self.barr * self.B) * coeff0**(-1) * coeff0_prime)

        theta1 = theta11 + theta12

        thetatot = theta0 + self.epsilon * theta1
    
        return thetatot

    def vr(self):
        coeff0 = self.C0(self.K)
        coeff0_prime = self.C0_prime(self.K)
        c3 = self.C3(self.K)
        

        _ , omega0 , _ = self.omega()

        vr0 = 0
        vr11 = omega0 * c3 * np.cos(omega0 * self.time)
        vr12 = ((1 - self.B) / (1 - self.B * self.barr)) * coeff0 * (self.barr**2 * self.B * coeff0 / (3 * (1 - self.B)) + 2 * coeff0_prime)
        vr1 = vr11 + vr12

        vrtot = vr0 + self.epsilon * vr1

        return vrtot    

if __name__== "__main__":
    par = dust_properties("silicate" , "CME" , "all" , "large")
    res = np.load("Files/rk45_t6_large_silicate_CMEsw.npz")
    x , y , _ , _ , _ , b , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    p = perturbed_functions(par , t , b , find_k = False)
    
    c0 = p.C0(p.K)
    k1 = p.K_calc()
    print(p.K)
    #print(f"r:" ,p.rad() , "beta:" ,beta ,"c0:", c0)
    #np.savez("Files/BKTEST_kbcomb1.npz" , b = b_arr , k = k_arr , r0 = r0_arr) #, betaval = betaval)
    #print(r0_arr[-1 , :])
    # plt.plot(t , k1)
    # plt.show()

    # ts = np.load("Files/BKTEST_kbcomb.npz")
    # b_arr , k_arr , r0_arr = [ts[i] for i in ("b" , "k" , "r0")]
    # print(k_arr[-1])
    

    # Bcst , Kcst , r0 = dust_cloud("Files/BKTEST_kbcomb1.npz")
    # print(Bcst[-1] , len(Kcst[-1 , :]) , len(r0[-1 , :]))
    
    # omegatot , rtot , theta , vr , c0 , beta = p.omega() , p.rad() , p.theta() , p.vr() , c0 , p.barr
    # omega , _ , _ = omegatot
    # r , _ , _ = rtot
    # plt.plot(t , np.sqrt(x**2+y**2) , label = "num")
    # plt.plot(t , r , label = "pert")
    # plt.legend()
    # plt.show()
    # np.savez("Files/multiscale_t5_079micron_silicate_slowsw.npz" , omega = omega , r = r , theta = theta , vr = vr , c0 = c0 , k = k1 , b = beta , t = p.time)
    #np.savez("Files/kvals_02micron_slowsw.npz" , k = p.K_calc() , t = t)

    