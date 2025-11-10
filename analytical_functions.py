import matplotlib.pyplot as plt
from constants import *
import numpy as np

"""calculates up to first order of beta hat"""
def betahat(t):
    """input: t (float), time to calculate for

       returns: betahat_sum (float), sum of all order terms """
       
    order_0 = betahat_0 
    order_1 = epsilon * t * betahat_0**2 / 3
    
    betahat_sum = order_0 + order_1
        
    return betahat_sum

"""acceleration formula, up to first order"""
def acceleration(x , y , t):
    """input: x (float), x position
              y (float), y position
              
       returns: ax , ay (tuple), acceleration in x and y direction"""
    
    r = np.sqrt(x**2 + y**2) 
    bhat = betahat(t)
    
    fraction = (1 - bhat * beta_0) / (1 - beta_0)
    
    ax = -fraction * x / r**3 #acc x
    ay = -fraction * y / r**3 #acc y

    return ax , ay    

    
    
    
    