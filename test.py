from scipy.stats import linregress
import numpy as np

def find_slope(x1 , x2 , y1 , y2):
    x = np.array([x1 , x2])
    y = np.array([y1 , y2])

    return linregress(x , y)

#Numerically lifetimes in years, changed eps units
true_lifetime_units = {"particle1": {
                 "size": 1.54079 * 10**(-6) , 
                 "silicate": {
                     "pr": {
                         "CME": 0.0 , #impact Sun
                         "slow": 0.0 , #impact Sun
                         "fast": 0.0 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #destroyed
                         "slow": 0.0 , #destroyed (linear approx)
                         "fast": 0.0 , #destroyed (linear approx)
                     } ,
                     "both": {
                         "CME": 0.0 , #destroyed
                         "slow": 0.0 , #impact Sun
                         "fast": 0.0 #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 0.0 , #impact Sun
                         "slow": 0.0 , #impact Sun
                         "fast": 0.0 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #destroyed (linear approx)
                         "slow": 0.0 , #destroyed (linear approx)
                         "fast": 0.0 #destroyed (linear approx)
                     } ,
                     "both": {
                         "CME": 0.0 , #impact Sun
                         "slow": 0.0 , #impact Sun
                         "fast": 0.0 #impact Sun
                     }
                 }} , 
                 
                 "particle2": {
                     "size": 0.01220 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 0.0 , #impact Sun
                         "slow": 0.0 , #impact Sun
                         "fast": 0.0 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #destroyed
                         "slow": 0.0 , #destroyed
                         "fast": 0.0 , #destroyed
                     } ,
                     "both": {
                         "CME": 0.0 , #destroyed
                         "slow": 0.0 , #destroyed
                         "fast": 0.0 #impact Sun
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
                         "CME": 0.0 , #impact Sun
                         "slow": 0.0 , #impact Sun
                         "fast": 0.0 #impact Sun
                     }
                 }}  , 
                 "particle3": {
                     "size": 0.00708 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 0.0 , #impact Sun
                         "slow": 0.0 , #impact Sun
                         "fast": 0.0 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #destroyed
                         "slow": 0.0 , #destroyed
                         "fast": 0.0 , #destroyed
                     } ,
                     "both": {
                         "CME": 0.0 , #destroyed
                         "slow": 0.0 , #destroyed
                         "fast": 0.0 #destroyed
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
                 }} }

#Slope calcs
x1 = true_lifetime_units.get("particle2" , {}).get("size") 
x2 = true_lifetime_units.get("particle3" , {}).get("size") 
x3 = true_lifetime_units.get("particle1" , {}).get("size") 

par2_sil_spu_slow = true_lifetime_units.get("particle2" , {}).get("silicate" , {}).get("sputtering" , {}).get("slow")
par2_sil_spu_fast = true_lifetime_units.get("particle2" , {}).get("silicate" , {}).get("sputtering" , {}).get("fast")

par2_car_spu_CME = true_lifetime_units.get("particle2" , {}).get("carbon" , {}).get("sputtering" , {}).get("CME")
par2_car_spu_slow = true_lifetime_units.get("particle2" , {}).get("carbon" , {}).get("sputtering" , {}).get("slow")
par2_car_spu_fast = true_lifetime_units.get("particle2" , {}).get("carbon" , {}).get("sputtering" , {}).get("fast")

par3_sil_spu_slow = true_lifetime_units.get("particle3" , {}).get("silicate" , {}).get("sputtering" , {}).get("slow")
par3_sil_spu_fast = true_lifetime_units.get("particle3" , {}).get("silicate" , {}).get("sputtering" , {}).get("fast")

par3_car_spu_CME = true_lifetime_units.get("particle3" , {}).get("carbon" , {}).get("sputtering" , {}).get("CME")
par3_car_spu_slow = true_lifetime_units.get("particle3" , {}).get("carbon" , {}).get("sputtering" , {}).get("slow")
par3_car_spu_fast = true_lifetime_units.get("particle3" , {}).get("carbon" , {}).get("sputtering" , {}).get("fast")

# sil_spu_slow_slope , sil_spu_slow_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_sil_spu_slow , par3_sil_spu_slow)
# par1_sil_spu_slow = sil_spu_slow_slope * x3 + sil_spu_slow_intercept
    
# sil_spu_fast_slope , sil_spu_fast_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_sil_spu_fast , par3_sil_spu_fast)
# par1_sil_spu_fast = sil_spu_fast_slope * x3 + sil_spu_fast_intercept

car_spu_CME_slope , car_spu_CME_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_car_spu_CME , par3_car_spu_CME)
par1_car_spu_CME = car_spu_CME_slope * x3 + car_spu_CME_intercept

car_spu_slow_slope , car_spu_slow_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_car_spu_slow , par3_car_spu_slow)
par1_car_spu_slow = car_spu_slow_slope * x3 + car_spu_slow_intercept

car_spu_fast_slope , car_spu_fast_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_car_spu_fast , par3_car_spu_fast)
par1_car_spu_fast = car_spu_fast_slope * x3 + car_spu_fast_intercept

if __name__ == "__main__":
    1
    # print(par1_car_spu_fast)
    

    