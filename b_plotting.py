import matplotlib.pyplot as plt
from constants import *
import numpy as np

"""function defining first order approximation of b for selected eps"""
def b(epsilon , t):
    """input: epsilon (float), cst value causing perturbation
              t (array), time vals in yr to evaluate b over
        
       returns: b (float), b as function of time for selected eps"""
       
    b0 = 1 - beta_0 #defining b0
    
    b = b0 + epsilon / 3 * (1 - b0)**2 * t #b expression
    
    return b

"""zeroth order x, y calcs"""
def xy_zerothorder(t):
    """input: t (array), time values used to calculate x0
        
       returns: pos (array), x0 zeroth order vals"""
       
    b0 = 1 - beta_0
    x0_hat = 1
    x0_abs = np.abs(x0_hat)
    
    alpha = b0**(1 / 2) * x0_abs**(-3 / 2)
    T = 1
    X = r0
    
    """0th order approx for x0 and y0"""
    
    x0 = X * np.cos(alpha * t) + X / T * np.sin(alpha * t)
    y0 = X * np.cos(alpha * t + np.pi / 2) + X / T * np.sin(alpha * t + np.pi / 2)
    
    pos = np.sqrt(x0**2 + y0**2) #find position to verify circular orbit
    
    return x0 , y0

"""first order approx to x and y"""
def xy_firstorder(x , y , t):
    x , y = xy_zerothorder(t)
    r = np.sqrt(x**2 + y**2)
    
    b0 = 1 - beta_0
    
    ax = - epsilon / 3 * (1 - b0) * t * x1 / r**3 
    
    
    
    
"""plotting b over selected time interval"""
def plotting(b , t_max):
    """input b (float), calculated b vals for each time
             t_max (float) , maximum time to calculate b"""
             
    plt.plot(t_max , b , label = "x0(t)")
    
    plt.xlabel("Time (yr)")
    plt.ylabel("x0(t)")
    plt.title("b as function of time, zeroth order approximation")
    plt.legend()
    plt.show()
    
if __name__ == "__main__":
    t_max = np.linspace(0 , 10 , 11)
    
    b = b(0.4 , t_max)
    
    
    