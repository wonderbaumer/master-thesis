import numpy as np
from config import eps , B , rhat0

"""calculates up to first order of beta hat from analytical solution"""
def betahat_analytical(t):
    """input: t (array), scaled time to calculate for

       returns: bvals (array), contains all orders as well as their sum"""
    
    frac = eps() / 3
    bvals = 1 / (1 - frac * t) #total expression
        
    return bvals

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




    