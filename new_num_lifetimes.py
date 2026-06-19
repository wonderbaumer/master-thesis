from scipy.stats import linregress
import numpy as np
from config import yr

"""Dictionaries calculating lifetime in years, based on constant sputtering and radially varying 
sputtering. Calculations have been done with a constant and varying sputtering, for a chosen set of
initial particle sizes. Simulations have been done for sputtering and PR effects combined as well as
separately by switching off the neglected effect during simulations. The resulting time has been 
multiplied by T to get physical time, and divided by yr to get the physical lifetime in years.

Some sputtering lifetimes are so long that the particles are not destroyed in simulated time,
for simplicity we use results from Baumann et al., 2020, where it is stated a liner relation between
all sputtering lifetimes for one solar wind condition and infer the longest sputtering lifetimes
using a linear increase from the smaller values along the same solar wind condition."""

"""Finding linear slope and intercept using two coordinates"""
def find_slope(x1 , x2 , y1 , y2):
    """input: x1 (float), x value 1
              x2 (float), x value 2
              y1 (float), y value 1
              y2 (float), y value 2

       returns: linregress(x , y) (tuple), slope, intercept and other params from provided coordinates"""
    
    x = np.array([x1 , x2])
    y = np.array([y1 , y2])

    return linregress(x , y)

"""Numerical lifetimes in years, constant sputtering"""
true_lifetime = {"C": {
                 "size": 1.54079 * 10**(-6) , 
                 "silicate": {
                     "pr": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "sputtering": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "both": {
                         "CME": 2.36730114e+03 , #destroyed
                         "slow": 3.20389148e+03 , #impact Sun
                         "fast": 3.22275984e+03 , #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "sputtering": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "both": {
                         "CME": 1.49282414e+03 , #impact Sun
                         "slow": 1.51390561e+03 , #impact Sun
                         "fast": 1.51422347e+03 , #impact Sun
                     }
                 }} , 
                 
                 "D": {
                     "size": 5.09598 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "sputtering": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "both": {
                         "CME": 6.16820043e+03 , #impact Sun
                         "slow": 1.10095397e+04 , #impact Sun
                         "fast": 1.10929794e+04 , #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 0.0 , 
                         "slow": 0.0 , 
                         "fast": 0.0 
                     } ,
                     "sputtering": {
                         "CME": 0.0 , 
                         "slow": 0.0 ,
                         "fast": 0.0 
                     } ,
                     "both": {
                         "CME": 5.49937775e+03 ,  #impact Sun
                         "slow": 5.78521536e+03 , #impact Sun
                         "fast": 5.78907350e+03 #impact Sun 
                     }
                 }} , 

                #  "C": {
                #      "size": 0.10165 * 10**(-6) ,
                #  "silicate": {
                #      "pr": {
                #          "CME": 0.0 ,
                #          "slow": 0.0 ,
                #          "fast": 0.0 ,
                #      } ,
                #      "sputtering": {
                #          "CME": 0.0 ,
                #          "slow": 0.0 ,
                #          "fast": 0.0 ,
                #      } ,
                #      "both": {
                #          "CME": 1.54741140e+02 , #destroyed
                #          "slow": 5.28990466e+02 , #impact Sun
                #          "fast": 5.54002031e+02 , #impacted Sun
                #      }
                #  } ,
                #  "carbon": {
                #      "pr": {
                #          "CME": 0.0 , 
                #          "slow": 0.0 , 
                #          "fast": 0.0 
                #      } ,
                #      "sputtering": {
                #          "CME": 0.0 , 
                #          "slow": 0.0 ,
                #          "fast": 0.0 
                #      } ,
                #      "both": {
                #          "CME": 0.0 , 
                #          "slow": 0.0 , 
                #          "fast": 0.0 
                #      }
                #  }} ,
                 "A": {
                     "size": 0.01220 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "sputtering": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "both": {
                         "CME": 1.72190836e+01 , #destroyed
                         "slow": 7.39183550e+02 , #destroyed
                         "fast": 1.38096192e+03 , #destroyed
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 0.0 , 
                         "slow": 0.0 , 
                         "fast": 0.0 
                     } ,
                     "sputtering": {
                         "CME": 0.0 , 
                         "slow": 0.0 ,
                         "fast": 0.0 
                     } ,
                     "both": {
                         "CME": 1.05418029e+02 , #impact Sun
                         "slow": 3.41296827e+02 , #impact Sun
                         "fast": 3.77842219e+02 #impact Sun 
                     }
                 }}}

"""Numerical lifetimes in years, variable sputtering"""
true_lifetime_variableeps = {"C": {
                 "size": 1.54079 * 10**(-6) , 
                 "silicate": {
                     "pr": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "sputtering": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "both": {
                         "CME": 3.26075235e+03 , #destroyed
                         "slow": 3.13064495e+03 , #impact Sun
                         "fast": 3.16541123e+03 , #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "sputtering": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "both": {
                         "CME": 0.0 , #yeet
                         "slow": 1.50585504e+03 , #impact Sun
                         "fast": 1.50646806e+03 , #impact Sun
                     }
                 }} , 
                 
                 "D": {
                     "size": 5.09598 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "sputtering": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "both": {
                         "CME": 5.86119732e+03 , #destroyed
                         "slow": 1.10095397e+04 , #impact Sun
                         "fast": 1.09967871e+04 , #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 0.0 , 
                         "slow": 0.0 , 
                         "fast": 0.0 
                     } ,
                     "sputtering": {
                         "CME": 0.0 , 
                         "slow": 0.0 ,
                         "fast": 0.0 
                     } ,
                     "both": {
                         "CME": 5.25520286e+03 , #impact Sun
                         "slow": 5.77881690e+03 , #impact Sun
                         "fast": 5.78656870e+03 , #impact SUn
                     }
                 }} , 

                #  "C": {
                #      "size": 0.10165 * 10**(-6) ,
                #  "silicate": {
                #      "pr": {
                #          "CME": 0.0 ,
                #          "slow": 0.0 ,
                #          "fast": 0.0 ,
                #      } ,
                #      "sputtering": {
                #          "CME": 0.0 ,
                #          "slow": 0.0 ,
                #          "fast": 0.0 ,
                #      } ,
                #      "both": {
                #          "CME": 4.49887270e+01 , #destroyed
                #          "slow": 4.77535600e+02 , #destroyed
                #          "fast": 5.21872800e+02 , #impact Sun
                #      }
                #  } ,
                #  "carbon": {
                #      "pr": {
                #          "CME": 0.0 , 
                #          "slow": 0.0 , 
                #          "fast": 0.0 
                #      } ,
                #      "sputtering": {
                #          "CME": 0.0 , 
                #          "slow": 0.0 ,
                #          "fast": 0.0 
                #      } ,
                #      "both": {
                #          "CME": 0.0 , 
                #          "slow": 0.0 , 
                #          "fast": 0.0 
                #      }
                #  }} ,
                 "A": {
                     "size": 0.01220 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "sputtering": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     } ,
                     "both": {
                         "CME": 1.70754274e+01 , #destroyed
                         "slow": 6.81124395e+02 , #destroyed
                         "fast": 1.19395152e+03 , #destroyed
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 0.0 , 
                         "slow": 0.0 , 
                         "fast": 0.0 
                     } ,
                     "sputtering": {
                         "CME": 0.0 , 
                         "slow": 0.0 ,
                         "fast": 0.0 
                     } ,
                     "both": {
                         "CME": 3.14274370e+01 , #destroyed
                         "slow": 2.96560310e+02 , #impact Sun
                         "fast": 3.51843168e+02 , #impact Sun
                     }
                 }}}