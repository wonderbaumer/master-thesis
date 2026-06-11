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
true_lifetime = {"A": {
                 "size": 1.54079 * 10**(-6) , 
                 "silicate": {
                     "pr": {
                         "CME": (1.9075403764698192e+04 * 0.17011521908343677) , #impact Sun
                         "slow": (1.9075403764698192e+04 * 0.17011521908343677) , #impact Sun
                         "fast": (1.9075403764698192e+04 * 0.17011521908343677) #impact Sun
                     } ,
                     "sputtering": {
                         "CME": (1.391586920435384e+04 * 0.17011521908343677) , #destroyed
                         "slow": (1.4512983980552342e+04 * 0.17011521908343677) , #assumed destroyed
                         "fast": (2.711353383139966e+04 * 0.17011521908343677) , #assumed destroyed
                     } ,
                     "both": {
                         "CME": (1.3915869204353734e+04 * 0.17011521908343677) , #destroyed
                         "slow": (1.883547822221666e+04 * 0.17011521908343677) , #impact Sun
                         "fast": (1.8946423881851584e+04 * 0.17011521908343677) #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": (8.155195188611925e+03 * 0.18571941895449673) , #impact Sun
                         "slow": (8.155195188611925e+03 * 0.18571941895449673) , #impact Sun
                         "fast": (8.155195188611925e+03 * 0.18571941895449673) #impact Sun
                     } ,
                     "sputtering": {
                         "CME": (1.5376640831096107e+04 * 0.18571941895449673) , #assumed destroyed
                         "slow": (7.333368763147648e+05 * 0.18571941895449673) , #assumed destroyed
                         "fast": (1.8751709349382098e+06 * 0.18571941895449673) #assumed destroyed
                     } ,
                     "both": {
                         "CME": (8.038775762306178e+03 * 0.18571941895449673) , #impact Sun
                         "slow": (8.152379519559163e+03 * 0.18571941895449673) , #impact Sun
                         "fast": (8.154092201989478e+03 * 0.18571941895449673) #impact Sun
                     }
                 }} , 
                 
                 "G": {
                     "size": 0.01220 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": (28958.819396637136e+04 * 0.16633430454094575) , #impact Sun
                         "slow": (28958.819396637136e+04 * 0.16633430454094575) , #impact Sun
                         "fast": (28958.819396637136e+04 * 0.16633430454094575) #impact Sun
                     } ,
                     "sputtering": {
                         "CME": (1.0352094029330188e+02 * 0.16633430454094575) , #destroyed
                         "slow": (4.443963329734232e+03 * 0.16633430454094575) , #destroyed
                         "fast": (8.30232778095182e+03 * 0.16633430454094575) #destroyed
                     } ,
                     "both": {
                         "CME": (1.0352094029330185e+02 * 0.16633430454094575) , #destroyed
                         "slow": (4.44396332973424e+03 * 0.16633430454094575) , #destroyed
                         "fast": (8.302327780951786e+03 * 0.16633430454094575) #destroyed
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": (401.1585920969129e+02 * 1.023791326777151) , #impact Sun
                         "slow": (401.1585920969129e+02 * 1.023791326777151) , #impact Sun
                         "fast": (401.1585920969129e+02 * 1.023791326777151) #impact Sun
                     } ,
                     "sputtering": {
                         "CME": (1.0928173000855521e+02 * 1.023791326777151) , #destroyed
                         "slow": (5.211822491208914e+03 * 1.023791326777151) , #destroyed
                         "fast": (1.3326832959341404e+04 * 1.023791326777151) #destroyed
                     } ,
                     "both": {
                         "CME": (1.0302896064960659e+02 * 1.023791326777151) , #impact Sun
                         "slow": (3.334806801435187e+02 * 1.023791326777151) , #impact Sun
                         "fast": (3.6919088304662847e+02 * 1.023791326777151) #impact Sun
                     }
                 }} , 
                 "H": {
                     "size": 0.00708 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": (3.0209106707628238e+04 * 0.16603575206890253) , #impact Sun
                         "slow": (3.0209106707628238e+04 * 0.16603575206890253) , #impact Sun
                         "fast": (3.0209106707628238e+04 * 0.16603575206890253) #impact Sun
                     } ,
                     "sputtering": {
                         "CME": (5.629813106023842e+01 * 0.16603575206890253) , #destroyed
                         "slow": (2.4167750916425866e+03 * 0.16603575206890253) , #destroyed
                         "fast": (4.515082032609089e+03 * 0.16603575206890253) #destroyed
                     } ,
                     "both": {
                         "CME": (5.629813106023846e+01 * 0.16603575206890253) , #destroyed
                         "slow": (2.4167750916425853e+03 * 0.16603575206890253) , #destroyed
                         "fast": (4.515082032609065e+03 * 0.16603575206890253) #destroyed
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": (5.889998280495578e+02 * 0.7158409348024901) , #impact Sun
                         "slow": (5.889998280495578e+02 * 0.7158409348024901) , #impact Sun
                         "fast": (5.889998280495578e+02 * 0.7158409348024901) #impact Sun
                     } ,
                     "sputtering": {
                         "CME": (8.484534787597285e+01 * 0.7158409348024901) , #destroyed
                         "slow": (4.0464118961131667e+03 * 0.7158409348024901) , #destroyed
                         "fast": (1.034683270874086e+04 * 0.7158409348024901) #destroyed
                     } ,
                     "both": {
                         "CME": (8.484534787597295e+01 * 0.7158409348024901) , #destroyed
                         "slow": (5.596875912887647e+02 * 0.7158409348024901) , #impact Sun
                         "fast": (5.761851971142123e+02 * 0.7158409348024901) #impact Sun
                     }
                 }} ,
                 "D": {
                     "size": 0.10165 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": (2.086304748095878e+03 * 0.28202776737574825) , #impact Sun
                         "slow": (2.086304748095878e+03 * 0.28202776737574825) , #impact Sun
                         "fast": (2.086304748095878e+03 * 0.28202776737574825) #impact Sun
                     } ,
                     "sputtering": {
                         "CME": (5.486734203335685e+02 * 0.28202776737574825) , #destroyed
                         "slow": (2.3553539535614473e+04 * 0.28202776737574825) , #destroyed
                         "fast": (4.400333466252357e+04 * 0.28202776737574825) #destroyed
                     } ,
                     "both": {
                         "CME": (5.486734203335675e+02 * 0.28202776737574825) , #destroyed
                         "slow": (1.8781639154328848e+03 * 0.28202776737574825) , #impact Sun
                         "fast": (1.9670131299476825e+03 * 0.28202776737574825) #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 0.0 , #impact Sun
                         "slow": 0.0 , #impact Sun
                         "fast": 0.0 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #destroyed
                         "slow": 0.0 , #destroyed
                         "fast": 0.0 #destroyed
                     } ,
                     "both": {
                         "CME": 0.0 , #destroyed
                         "slow": 0.0 , #impact Sun
                         "fast": 0.0 #impact Sun
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
true_lifetime_variableeps = {"A": {
                 "size": 1.54079 * 10**(-6) , 
                 "silicate": {
                     "pr": {
                         "CME": 0.0 ,
                         "slow": 0.0 ,
                         "fast": 0.0 
                     } ,
                     "sputtering": {
                         "CME": 0.0 , 
                         "slow": 0.0 , 
                         "fast": 0.0 , 
                     } ,
                     "both": {
                         "CME": 3.2799220179287167e+03 , #destroyed 
                         "slow": 3.165877728790449e+03 , #impact Sun 
                         "fast": 3.2018350461736945e+03 #impact Sun
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
                         "CME": 0.0 , #yeet
                         "slow": 1.513539876418186e+03 , #impact Sun
                         "fast": 1.5141695995692314e+03 #impact Sun
                     }
                 }} , 
                 
                 "G": {
                     "size": 0.01220 * 10**(-6) ,
                 "silicate": {
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
                         "CME": 1.7076197230341602e+01 , #destroyed
                         "slow": 6.822441704258237e+02 , #destroyed
                         "fast": 1.197496850939118e+03 #destroyed
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
                         "CME": 3.1549098696479202e+01 , #destroyed
                         "slow": 2.966635951570946e+02 , #impact Sun
                         "fast": 3.519646536759617e+02 #impact Sun
                     }
                 }} , 
                 "H": {
                     "size": 0.00708 * 10**(-6) ,
                 "silicate": {
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
                         "CME": 9.324368818449967 , #destroyed 
                         "slow": 3.8514844692983246e+02 , #destroyed
                         "fast": 6.956122480723823e+02 #destroyed
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
                         "CME": 4.4403387296926226e+01 , #destroyed 
                         "slow": 3.846379010359448e+02 , #destroyed
                         "fast": 4.0476340615297164e+02 #impact Sun
                     }
                 }} ,
                 "D": {
                     "size": 0.10165 * 10**(-6) ,
                 "silicate": {
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
                         "CME": 4.4995574900932155e+01 , #destroyed 
                         "slow": 4.781054446295344e+02 , #destroyed
                         "fast": 5.226594334458505e+02 #destroyed
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
                 }} }

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

    