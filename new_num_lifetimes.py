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
true_lifetime = {"D": {
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
                 
                 "E": {
                     "size": 40.22706 * 10**(-6) ,
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
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
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
                         "CME": 0.0 , 
                         "slow": 0.0 , 
                         "fast": 0.0 
                     }
                 }} , 

                 "C": {
                     "size": 0.10165 * 10**(-6) ,
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
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
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
                         "CME": 0.0 , 
                         "slow": 0.0 , 
                         "fast": 0.0 
                     }
                 }} ,
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
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
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
                         "CME": 0.0 , 
                         "slow": 0.0 , 
                         "fast": 0.0 
                     }
                 }}}

"""Slope and intercept calcs, constant sputtering. Using y = ax + b"""
# x1 = true_lifetime.get("G" , {}).get("size") 
# x2 = true_lifetime.get("H" , {}).get("size") 
# x3 = true_lifetime.get("A" , {}).get("size") 

# G_sil_spu_slow = true_lifetime.get("G" , {}).get("silicate" , {}).get("sputtering" , {}).get("slow")
# G_sil_spu_fast = true_lifetime.get("G" , {}).get("silicate" , {}).get("sputtering" , {}).get("fast")

# G_car_spu_CME = true_lifetime.get("G" , {}).get("carbon" , {}).get("sputtering" , {}).get("CME")
# G_car_spu_slow = true_lifetime.get("G" , {}).get("carbon" , {}).get("sputtering" , {}).get("slow")
# G_car_spu_fast = true_lifetime.get("G" , {}).get("carbon" , {}).get("sputtering" , {}).get("fast")

# H_sil_spu_slow = true_lifetime.get("H" , {}).get("silicate" , {}).get("sputtering" , {}).get("slow")
# H_sil_spu_fast = true_lifetime.get("H" , {}).get("silicate" , {}).get("sputtering" , {}).get("fast")

# H_car_spu_CME = true_lifetime.get("H" , {}).get("carbon" , {}).get("sputtering" , {}).get("CME")
# H_car_spu_slow = true_lifetime.get("H" , {}).get("carbon" , {}).get("sputtering" , {}).get("slow")
# H_car_spu_fast = true_lifetime.get("H" , {}).get("carbon" , {}).get("sputtering" , {}).get("fast")

# sil_spu_slow_slope , sil_spu_slow_intercept , _ , _ , _ = find_slope(x1 , x2 , G_sil_spu_slow 
#                                                                      , H_sil_spu_slow)
# A_sil_spu_slow = sil_spu_slow_slope * x3 + sil_spu_slow_intercept
    
# sil_spu_fast_slope , sil_spu_fast_intercept , _ , _ , _ = find_slope(x1 , x2 , G_sil_spu_fast 
#                                                                      , H_sil_spu_fast)
# A_sil_spu_fast = sil_spu_fast_slope * x3 + sil_spu_fast_intercept

# car_spu_CME_slope , car_spu_CME_intercept , _ , _ , _ = find_slope(x1 , x2 , G_car_spu_CME 
#                                                                    , H_car_spu_CME)
# A_car_spu_CME = car_spu_CME_slope * x3 + car_spu_CME_intercept

# car_spu_slow_slope , car_spu_slow_intercept , _ , _ , _ = find_slope(x1 , x2 , G_car_spu_slow 
#                                                                      , H_car_spu_slow)
# A_car_spu_slow = car_spu_slow_slope * x3 + car_spu_slow_intercept

# car_spu_fast_slope , car_spu_fast_intercept , _ , _ , _ = find_slope(x1 , x2 , G_car_spu_fast 
#                                                                      , H_car_spu_fast)
# A_car_spu_fast = car_spu_fast_slope * x3 + car_spu_fast_intercept

"""Numerical lifetimes in years, variable sputtering"""
true_lifetime_variableeps = {"D": {
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
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
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
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
                     }
                 }} , 
                 
                 "E": {
                     "size": 50 * 10**(-6) ,
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
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
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
                         "CME": 0.0 , 
                         "slow": 0.0 , 
                         "fast": 0.0 
                     }
                 }} , 

                 "C": {
                     "size": 0.10165 * 10**(-6) ,
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
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
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
                         "CME": 0.0 , 
                         "slow": 0.0 , 
                         "fast": 0.0 
                     }
                 }} ,
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
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 ,
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
                         "CME": 0.0 , 
                         "slow": 0.0 , 
                         "fast": 0.0 
                     }
                 }}}

"""Slope and intercept calcs, variable sputtering. Using y = ax + b"""
# x11 = true_lifetime.get("E" , {}).get("size") 
# x21 = true_lifetime.get("F" , {}).get("size") 
# x31 = true_lifetime.get("A" , {}).get("size") 

# par2_sil_spu_slow1 = true_lifetime.get("E" , {}).get("silicate" , {}).get("sputtering" , {}).get("slow")
# par2_sil_spu_fast1 = true_lifetime.get("E" , {}).get("silicate" , {}).get("sputtering" , {}).get("fast")

# par2_car_spu_CME1 = true_lifetime.get("E" , {}).get("carbon" , {}).get("sputtering" , {}).get("CME")
# par2_car_spu_slow1 = true_lifetime.get("E" , {}).get("carbon" , {}).get("sputtering" , {}).get("slow")
# par2_car_spu_fast1 = true_lifetime.get("E" , {}).get("carbon" , {}).get("sputtering" , {}).get("fast")

# par3_sil_spu_slow1 = true_lifetime.get("F" , {}).get("silicate" , {}).get("sputtering" , {}).get("slow")
# par3_sil_spu_fast1 = true_lifetime.get("F" , {}).get("silicate" , {}).get("sputtering" , {}).get("fast")

# par3_car_spu_CME1 = true_lifetime.get("F" , {}).get("carbon" , {}).get("sputtering" , {}).get("CME")
# par3_car_spu_slow1 = true_lifetime.get("F" , {}).get("carbon" , {}).get("sputtering" , {}).get("slow")
# par3_car_spu_fast1 = true_lifetime.get("F" , {}).get("carbon" , {}).get("sputtering" , {}).get("fast")

# sil_spu_slow_slope1 , sil_spu_slow_intercept1 , _ , _ , _ = find_slope(x11 , x21 , par2_sil_spu_slow1 , par3_sil_spu_slow1)
# par1_sil_spu_slow1 = sil_spu_slow_slope1 * x31 + sil_spu_slow_intercept1
    
# sil_spu_fast_slope1 , sil_spu_fast_intercept1 , _ , _ , _ = find_slope(x11 , x21 , par2_sil_spu_fast1 , par3_sil_spu_fast1)
# par1_sil_spu_fast1 = sil_spu_fast_slope1 * x31 + sil_spu_fast_intercept1

# car_spu_CME_slope1 , car_spu_CME_intercept1 , _ , _ , _ = find_slope(x11 , x21 , par2_car_spu_CME1 , par3_car_spu_CME1)
# par1_car_spu_CME1 = car_spu_CME_slope1 * x31 + car_spu_CME_intercept1

# car_spu_slow_slope1 , car_spu_slow_intercept1 , _ , _ , _ = find_slope(x11 , x21 , par2_car_spu_slow1 , par3_car_spu_slow1)
# par1_car_spu_slow1 = car_spu_slow_slope1 * x31 + car_spu_slow_intercept1

# car_spu_fast_slope1 , car_spu_fast_intercept1 , _ , _ , _ = find_slope(x11 , x21 , par2_car_spu_fast1 , par3_car_spu_fast1)
# par1_car_spu_fast1 = car_spu_fast_slope1 * x31 + car_spu_fast_intercept1