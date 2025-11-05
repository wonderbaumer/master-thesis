import matplotlib.pyplot as plt
from constants import *
import numpy as np

"""calculates up to first order of beta hat"""
def betahat(t , eps):
    """input: t (float), time to calculate for
              eps (float), perturbation term
        
       returns: betahat_sum (float), sum of all order terms """
       
    order_0 =  1 / betahat_0 
    order_1 = eps * t / (3 * betahat_0**2)
    
    betahat_sum = order_0 + order_1
    
    return betahat_sum

   
    
"""plotting b over selected time interval"""
def plotting(y , t_max):
    """input b (float), calculated b vals for each time
             t_max (float) , maximum time to calculate b"""
             
    plt.plot(t_max , y , label = "betahat(t), eps=0.019")
    
    plt.xlabel("Time (yr)")
    plt.ylabel("x0(t)")
    plt.title("beta hat as function of time")
    plt.legend()
    plt.show()
    
if __name__ == "__main__":
    t0 = 0 #initial time in s
    #t_tot = 3.16e10 #total time in s
    
    t_max = np.linspace(t0 , 1 , 100)
    epsilon = 0.019
    
    y = betahat(t_max , epsilon)
    
    plotting(y , t_max)
    
    
    