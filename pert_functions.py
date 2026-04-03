import numpy as np
from config import sil_beta , car_beta , init_vals 
from scipy.constants import c

class perturbed_functions():

    def __init__(self , particle , t):

        self.particle = particle
        self.material = particle.material
        self.size = particle.size

        self.r = particle.r
        self.B = particle.B
        self.V = particle.V
        
        #dt , t_tot = t
        #self.time = np.arange(0 , t_tot , dt)
        self.time = t
        self.epsilon = particle.eps()
        self.K = self.V / c / self.epsilon
        self.coeff3 = self.C3()
        self.d0 = self.D0()

        

    """calculates up to first order of beta hat from analytical solution"""
    def betahat_analytical(self):
        """input: t (array), scaled time to calculate for

        returns: bvals (array), contains all orders as well as their sum"""
    
        bvals = 1 / (1 -  self.epsilon * self.time / 3) #total expression
        
        return bvals

    ###WITH DRAG###
    def C0(self , b_func):

        first = -4 * self.B * self.K / (1 - self.B)**3
        second = (3 * b_func**4 * self.B / 4 - 4 * self.B**3 * b_func**3 + 9 * self.B**2 * b_func**2 - 12 * self.B * b_func + 3 * np.log(b_func))
        third = 1 + (4 * self.B * self.K / (1 - self.B)**3) * (9 * self.B**2 - 45 * self.B / 4 - 4 * self.B**3)

        tot = (first * second + third)**(1 / 4)

        return tot
    
    def C0_prime(self , b_func):

        numerator = -self.B * self.K * (-4 * self.B**3 * b_func**4 + 6 * self.B**2 * b_func**3 - 4 * self.B * b_func**2 + self.B * b_func**5 + b_func)
        denom1 = 4 * self.B * self.K * (-4 * self.B**3 + 9 * self.B**2 - 45 * self.B / 4)
        denom2 = -4 * self.B * self.K * (-4 * self.B**3 * b_func**3 + 9 * self.B**2 * b_func**2 - 12 * self.B * b_func + 3 * self.B * b_func**4 / 4 + 3 * np.log(b_func))
        denom3 = (1 - self.B)**(-3)
        tot_denom = (denom1 + denom2 + denom3)**(3 / 4)

        tot = numerator / tot_denom
        
        return tot
    
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
        coeff0_prime = self.C0_prime(b_func)

        _ , omega0 , _ = self.omega()
        _ , r0 , _ = self.rad()

        theta0 = ((1 - self.B) / (1 - self.B * b_func))**(-2) * coeff0**(-3) * self.time + self.d0
        theta11 = -2 / r0 * self.coeff3 * np.cos(omega0 * self.time)
        theta12 = -self.time**2 / 2 * coeff0**(-3) * ((1 - b_func * self.B) / (1 - self.B)**2) * (-2 * self.B * b_func**2 / 3 - 3 * (1 - b_func * self.B) * coeff0**(-3) * coeff0_prime)

        theta1 = theta11 + theta12

        thetatot = theta0 + self.epsilon * theta1
    
        return thetatot

    def vr(self):
        b_func = self.betahat_analytical()
        coeff0 = self.C0(b_func)
        coeff0_prime = self.C0_prime(b_func)

        _ , omega0 , _ = self.omega()

        vr0 = 0
        vr11 = omega0 * self.coeff3 * np.cos(omega0 * self.time)
        vr12 = ((1 - self.B) / (1 - self.B * b_func)) * coeff0 * (b_func**2 * self.B * coeff0 / (3 * (1 - self.B)) + 2 * coeff0_prime)
        vr1 = vr11 + vr12

        vrtot = vr0 + self.epsilon * vr1

        return vrtot


if __name__== "__main__":
    1
    #particle = dust_properties(material, sw, species, size)

    