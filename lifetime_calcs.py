from scipy.stats import linregress
import numpy as np

def find_slope(x1 , x2 , y1 , y2):
    x = np.array([x1 , x2])
    y = np.array([y1 , y2])

    return linregress(x , y)

#Numerically lifetimes in years
true_lifetime = {"particle1": {
                 "size": 1.54079 * 10**(-6) , 
                 "silicate": {
                     "pr": {
                         "CME": 3.245016490536648e+03 , #impact Sun
                         "slow": 3.245016490536648e+03 , #impact Sun
                         "fast": 3.245016490536648e+03 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 3.380760903403005e+02 , #destroyed
                         "slow": 1.4512983980552342e+04 , #assumed destroyed
                         "fast": 2.711353383139966e+04 , #assumed destroyed
                     } ,
                     "both": {
                         "CME": 3.380760903395063e+02 , #destroyed
                         "slow": 2.976750343468511e+03 , #impact Sun
                         "fast": 3.0965710813731766e+03 #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 1.514578111889514e+03 , #impact Sun
                         "slow": 1.514578111889514e+03 , #impact Sun
                         "fast": 1.514578111889514e+03 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 1.5376640831096107e+04 , #assumed destroyed
                         "slow": 7.333368763147648e+05 , #assumed destroyed
                         "fast": 1.8751709349382098e+06 #assumed destroyed
                     } ,
                     "both": {
                         "CME": 1.4929509282623908e+03 , #impact Sun
                         "slow": 1.5140550181222495e+03 , #impact Sun
                         "fast": 1.5143731994043846e+03 #impact Sun
                     }
                 }} , 
                 
                 "particle2": {
                     "size": 0.01220 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 4.816845084666488e+03 , #impact Sun
                         "slow": 4.816845084666488e+03 , #impact Sun
                         "fast": 4.816845084666488e+03 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 2.4590705303254685 , #destroyed
                         "slow": 1.0556336940893623e+02 , #destroyed
                         "fast": 1.9721623007791672e+02 #destroyed
                     } ,
                     "both": {
                         "CME": 2.4590705303267577 , #destroyed
                         "slow": 1.0556336940893613e+02 , #destroyed
                         "fast": 1.9721623007791501e+02 #destroyed
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 4.1070268725095247e+02 , #impact Sun
                         "slow": 4.1070268725095247e+02 , #impact Sun
                         "fast": 4.1070268725095247e+02 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 1.118453667761682e+02 , #destroyed
                         "slow": 5.334086475899605e+03 , #destroyed
                         "fast": 1.3639466726831513e+04 #destroyed
                     } ,
                     "both": {
                         "CME": 1.0546861796842475e+02 , #impact Sun
                         "slow": 3.413983230000514e+02 , #impact Sun
                         "fast": 3.7796525897239525e+02 #impact Sun
                     }
                 }} , 
                 "particle3": {
                     "size": 0.00708 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 5.015791751530783e+03 , #impact Sun
                         "slow": 5.015791751530783e+03 , #impact Sun
                         "fast": 5.015791751530783e+03 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 1.3349240021964774 , #destroyed
                         "slow": 5.730582910770815e+01 , #destroyed
                         "fast": 1.0706023918515456e+02 #destroyed
                     } ,
                     "both": {
                         "CME": 1.3349240021963173 , #destroyed
                         "slow": 5.73058291077082e+01 , #destroyed
                         "fast": 1.0706023918515443e+02 #destroyed
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 4.2163018750950135e+02 , #impact Sun
                         "slow": 4.2163018750950135e+02 , #impact Sun
                         "fast": 4.2163018750950135e+02 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 6.071605624991979e+01 , #destroyed
                         "slow": 2.8956469440598116e+03 , #destroyed
                         "fast": 7.4042819374228075e+03 #destroyed
                     } ,
                     "both": {
                         "CME": 6.071605624992002e+01 , #destroyed
                         "slow": 4.00641455738426e+02 , #impact Sun
                         "fast": 4.124541673429061e+02 #impact Sun
                     }
                 }} ,
                 "particle4": {
                     "size": 0.10165 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 5.883958701709033e+02 , #impact Sun
                         "slow": 5.883958701709033e+02 , #impact Sun
                         "fast": 5.883958701709033e+02 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 2.2098700792098164e+01 , #destroyed
                         "slow": 9.486565295544124e+02 , #destroyed
                         "fast": 1.7723047819055546e+03 #destroyed
                     } ,
                     "both": {
                         "CME": 2.20987007920982e+01 , #destroyed
                         "slow": 3.695179273496882e+02 , #impact Sun
                         "fast": 4.2996592938159836e+02 #impact Sun
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

#Slope calcs
x1 = true_lifetime.get("particle2" , {}).get("size") 
x2 = true_lifetime.get("particle3" , {}).get("size") 
x3 = true_lifetime.get("particle1" , {}).get("size") 

par2_sil_spu_slow = true_lifetime.get("particle2" , {}).get("silicate" , {}).get("sputtering" , {}).get("slow")
par2_sil_spu_fast = true_lifetime.get("particle2" , {}).get("silicate" , {}).get("sputtering" , {}).get("fast")

par2_car_spu_CME = true_lifetime.get("particle2" , {}).get("carbon" , {}).get("sputtering" , {}).get("CME")
par2_car_spu_slow = true_lifetime.get("particle2" , {}).get("carbon" , {}).get("sputtering" , {}).get("slow")
par2_car_spu_fast = true_lifetime.get("particle2" , {}).get("carbon" , {}).get("sputtering" , {}).get("fast")

par3_sil_spu_slow = true_lifetime.get("particle3" , {}).get("silicate" , {}).get("sputtering" , {}).get("slow")
par3_sil_spu_fast = true_lifetime.get("particle3" , {}).get("silicate" , {}).get("sputtering" , {}).get("fast")

par3_car_spu_CME = true_lifetime.get("particle3" , {}).get("carbon" , {}).get("sputtering" , {}).get("CME")
par3_car_spu_slow = true_lifetime.get("particle3" , {}).get("carbon" , {}).get("sputtering" , {}).get("slow")
par3_car_spu_fast = true_lifetime.get("particle3" , {}).get("carbon" , {}).get("sputtering" , {}).get("fast")

sil_spu_slow_slope , sil_spu_slow_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_sil_spu_slow , par3_sil_spu_slow)
par1_sil_spu_slow = sil_spu_slow_slope * x3 + sil_spu_slow_intercept
    
sil_spu_fast_slope , sil_spu_fast_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_sil_spu_fast , par3_sil_spu_fast)
par1_sil_spu_fast = sil_spu_fast_slope * x3 + sil_spu_fast_intercept

car_spu_CME_slope , car_spu_CME_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_car_spu_CME , par3_car_spu_CME)
par1_car_spu_CME = car_spu_CME_slope * x3 + car_spu_CME_intercept

car_spu_slow_slope , car_spu_slow_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_car_spu_slow , par3_car_spu_slow)
par1_car_spu_slow = car_spu_slow_slope * x3 + car_spu_slow_intercept

car_spu_fast_slope , car_spu_fast_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_car_spu_fast , par3_car_spu_fast)
par1_car_spu_fast = car_spu_fast_slope * x3 + car_spu_fast_intercept

if __name__ == "__main__":
    
    1
    