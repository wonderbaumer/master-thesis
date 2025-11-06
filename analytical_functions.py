import matplotlib.pyplot as plt
from constants import *
import numpy as np

"""calculates up to first order of beta hat"""
def betahat(t , eps):
    """input: t (float), time to calculate for
              eps (float), perturbation term
        
       returns: betahat_sum (float), sum of all order terms """
       
    order_0 = betahat_0 
    order_1 = eps * t * betahat_0**2 / 3
    
    betahat_sum = order_0 + order_1
    
    beta = np.array([order_0 , order_1 , betahat_sum])
    
    return beta

"""zeroth order radial acceleration"""
def radial_0(r , theta_vel , beta):
    """input: r (float), radial distance
              theta_vel (float), angular velocity
              beta (float), zeroth order beta value
       returns: r_0 (array), zeroth order approximated radial distance"""
       
    L = r **2 * theta_vel
    
    term1 = 2 / ((1 - beta[0]) * r**2)
    term2 = L**2 / r**3
    
    r_0 = np.array([term1 + term2])
    
    return r_0
   
    
"""plotting b over selected time interval"""
def plotting(y , t_max):
    """input b (float), calculated b vals for each time
             t_max (float) , maximum time to calculate b"""
             
    plt.plot(t_max / yr , y , label = "betahat(t), eps=0.019")
    
    plt.xlabel("Time (yr)")
    plt.ylabel("beta hat")
    plt.title("beta hat as function of time")
    plt.legend()
    plt.show()
    
if __name__ == "__main__":
    t0 = 0 #initial time in s
    t_tot = yr / yr #total time in s
    
    t_max = np.linspace(t0 , t_tot , 9978)
    epsilon = 0.019
    
    
    