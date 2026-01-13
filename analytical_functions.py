import numpy as np
from constants import *
from polar_to_cart import *
from config import *

"""calculates up to first order of beta hat from explicit solution"""
def betahat_analytical(t):
    """input: t (array), scaled time to calculate for

       returns: bvals (array), contains all orders as well as their sum"""
       
    order_0 = np.ones_like(t) #zeroth order
    order_1 = eps() * t / 3 #first order
    tot = order_0 + order_1 #both orders

    bvals = np.column_stack((order_0 , order_1 , tot)) 
        
    return bvals

"""calculates up to first order of betahat from perturbed expression"""
def betahat_pert(t):
    """input: t (array), scaled time to evaluate

       returns: bvals (array), contain all orders as well as sum of orders """
    
    order_0 = 1 * np.ones_like(t) #zeroth order
    order_1 = eps() * t / 3 #first order
    tot = order_0 + order_1 #sum of orders

    bvals = np.column_stack((order_0 , order_1 , tot))

    return bvals

"""calculates radial position up to first order in epsilon"""
def radial_position(t):
    """input: t (array), time (scaled) to calculate for

       returns: r (array), radial position up to first order in epsilon"""
    #print(-eps()*(beta0 * np.sin(t) / (3 * (1 - beta0))), eps()*( beta0 * t /(3 * (1 - beta0))))
    r0 = rhat0 #zeroth order
    r1 = -beta0 * np.sin(t) / (3 * (1 - beta0)) + beta0 * t /(3 * (1 - beta0)) #first order

    r = r0 + eps() * r1 #total expansion

    return r

"""calculates angular position up to first order in epsilon"""
def angular_position(t):
    """input: t (array), time (scaled) to calculate for

       returns: theta (array), angular position up to first order in epsilon"""
    #print(-beta0 * t**2 / (3 * (1 - beta0))*eps(), -eps()*2 * beta0 * np.cos(t) / (3 * (1 - beta0)),t)
    theta0 = t #zeroth order
    theta1 = -beta0 * t**2 / (3 * (1 - beta0))  - 2 * beta0 * np.cos(t) / (3 * (1 - beta0)) #first order
    theta = theta0 + eps() * theta1 #total expansion

    return theta

"""calculates the total, analytical orbit from r and theta perturbation expressions, both up to first
order in epsilon"""
def analytical_orbit(t):
    """input: t (array) , scaled time over which to plot the orbit
        
       returns: x,y (tuple), x and y coordinates"""

    r = radial_position(t) #r calcs
    theta = angular_position(t) #theta calcs
    x , y = r * np.cos(theta) , r * np.sin(theta) #converting to x and y

    return x , y
 
if __name__ == "__main__":
    dt , t_tot = t4
    that = np.arange(0 , t_tot , dt) / T

    #angular_position(that)
    radial_position(that)

    #print(angular_position(that))