import numpy as np
from config import sil_beta , car_beta , init_vals

class perturbed_functions():

    def __init__(self , material , size , drag_coeff , particle , t):
        self.material = material
        self.size = size
        self.r = init_vals[self.size]["r"]
        self.B = init_vals[self.size]["B"][self.material.lower()]
        self.K = drag_coeff
        self.particle = particle
        self.t = t
        self.eps = self.particle.eps()
        self.b_func = self.betahat_analytical()
        self.c = self.C0()

    """calculates up to first order of beta hat from analytical solution"""
    def betahat_analytical(self):
        """input: t (array), scaled time to calculate for

        returns: bvals (array), contains all orders as well as their sum"""
    
        bvals = 1 / (1 -  self.eps * self.t / 3) #total expression
        
        return bvals

    ###WITH DRAG###
    def C0(self):

        first = -4 * self.B * self.K / (1 - self.B)**3
        second = 3 * self.b_func**4 * self.B / 4 - 4 * self.B**3 * self.b_func**3 + 9 * self.B**2 * self.b_func**2 - 12 * self.B * self.b_func + 3 * np.log(self.b_func)
        third = 1 + 4 * self.B * self.K / (1 - self.B)**3 * (9 * self.B**2 - 45 * self.B / 4 - 4 * self.B**3)

        tot = first * second + third

        return tot**(1 / 4)

    def omega(self):
        omega0 = ((1 - self.B) / (1 - self.b_func * self.B))**(-2) * self.c**(-3)
        omega1 = -2 * omega0 / (((1 - self.B) / (1 - self.b_func * self.B)) * self.c**2) * (-2 * self.K - self.B / (3 * (1 - self.B))) * np.sin(omega0 * self.t)

        omegatot = omega0 + self.eps * omega1

        return omegatot , omega0 , omega1

def r(t , beta_func , c , om):
    _ , omega0 , _ = om

    r0 = ((1 - B) / (1 - beta_func * B)) * c**2
    r1 = (-2 * K - B / (3 * (1 - B))) * np.sin(omega0 * t)

    rtot = r0 + eps() * r1

    return rtot

def theta(t , beta_func):
    theta0 = t
    theta1 = -2 * B /(9 * (1 - B)) * (1 - ((1 - beta_func * B) / (1 - B))**(-1 / 3) * np.cos(t))

    thetatot = theta0 + eps() * theta1
    
    return thetatot

def vr(t , beta_func):
    vr0 = 0
    vr1 = B / (9 * (1 - B)) * (np.cos(t) - beta_func**2 * ((1 - beta_func * B) / (1 - B))**(-2 / 3))

    vrtot = vr0 + eps() * vr1

    return vrtot

###WITHOUT DRAG###

"""calculates up to first order of betahat from perturbed expression"""
def betahat_pert(t):
    """input: t (array), scaled time to evaluate

       returns: bvals (array), contain all orders as well as sum of orders """
    
    betahat_0 = 1 * np.ones_like(t) #zeroth order
    betahat_1 = eps() * t / 3 #first order
    betahat_tot = betahat_0 + betahat_1 #sum of orders

    return betahat_tot

"""calculates radial position up to first order in epsilon"""
def rhat_pert(t):
    """input: t (array), time (scaled) to calculate for

       returns: r (array), radial position up to first order in epsilon"""

    rhat_0 = rhat0 #zeroth order
    rhat_1 = -B * np.sin(t) / (3 * (1 - B)) + B * t /(3 * (1 - B)) #first order
    rhat = rhat_0 + eps() * rhat_1 #total expansion

    return rhat

"""calculates angular position up to first order in epsilon, perturbed expression"""
def thetahat_pert(t):
    """input: t (array), \time (scaled) to calculate for

       returns: theta (array), angular position up to first order in epsilon"""

    thetahat_0 = t #zeroth order 
    thetahat_1 = -B * t**2 / (3 * (1 - B))  - 2 * B * np.cos(t) / (3 * (1 - B)) #first order
    thetahat_tot = thetahat_0 + eps() * thetahat_1 #total expansion

    return thetahat_tot

"""calculates the total perturbed orbit in x and y coordinates"""
def perturbed_orbit(t):
    """input: t (array) , scaled time over which to plot the orbit
        
       returns: x,y (tuple), x and y coordinates"""

    rhat = rhat_pert(t) #r calcs
    thetahat = thetahat_pert(t) #theta calcs
    x , y = rhat * np.cos(thetahat) , rhat * np.sin(thetahat) #converting to x and y

    return x , y
 
"calculates perturbed velocity based on t hat"
def vrhat_pert(t):
    """input: t (array), t hat
    
       returns: v (array), v hat values"""
    
    vhat_0 = 0 #zeroth order v hat
    vhat_1 = B / (3 * (1 - B)) - B * np.cos(t) / (3 * (1 - B)) #first order v hat
    vhat_tot = vhat_0 + eps() * vhat_1 #tot v hat up to first order in epsilon

    return vhat_tot

"calculates perturbed angular velocity based on t hat"
def omegahat_pert(t):
    """input: t (array), t hat
       
       returns: omega (array), omega hat values"""
    
    omegahat_0 = 1 #zeroth order omega hat
    omegahat_1 = 2 * B * np.sin(t) / (3 * (1 - B)) - 2 * B * t / (3 * (1 - B)) #first order omega hat
    omegahat_tot = omegahat_0 + eps() * omegahat_1 #tot omega hat up to first order in epsilon

    return omegahat_tot

if __name__== "__main__":
    1
    #particle = dust_properties(material, sw, species, size)

    