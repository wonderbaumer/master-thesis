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
                         "CME": 0.0 , #destroyed
                         "slow": 0.0 , #destroyed (linear approx)
                         "fast": 0.0 , #destroyed (linear approx)
                     } ,
                     "both": {
                         "CME": 7.896128519889596e+02 , #destroyed
                         "slow": 3.203314127759356e+03 , #impact Sun
                         "fast": 3.222822151045436e+03 #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 1.514578111889514e+03 , #impact Sun
                         "slow": 1.514578111889514e+03 , #impact Sun
                         "fast": 1.514578111889514e+03 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #destroyed (linear approx)
                         "slow": 0.0 , #destroyed (linear approx)
                         "fast": 0.0 #destroyed (linear approx)
                     } ,
                     "both": {
                         "CME": 1.491727867903059e+03, #impact Sun
                         "slow": 1.5140544689126152e+03 , #impact Sun
                         "fast": 1.5143731557036722e+03 #impact Sun
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
                         "CME": 0.0 , #destroyed
                         "slow": 0.0 , #destroyed
                         "fast": 0.0 , #destroyed
                     } ,
                     "both": {
                         "CME": 6.248724151259423 , #destroyed
                         "slow": 2.6824622059212913e+02 , #destroyed
                         "fast": 5.0114456040989717e+02 #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 4.1070268725095247e+02 , #impact Sun
                         "slow": 4.1070268725095247e+02 , #impact Sun
                         "fast": 4.1070268725095247e+02 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #destroyed
                         "slow": 0.0 , #destroyed
                         "fast": 0.0 #destroyed
                     } ,
                     "both": {
                         "CME": 4.060133615399969e+01 , #impact Sun
                         "slow": 3.393155772855181e+02 , #impact Sun
                         "fast": 3.7748031012755035e+02 #impact Sun
                     }
                 }}  , 
                 "particle3": {
                     "size": 0.00708 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 5.015791751530783e+03 , #impact Sun
                         "slow": 5.015791751530783e+03 , #impact Sun
                         "fast": 5.015791751530783e+03 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #destroyed
                         "slow": 0.0 , #destroyed
                         "fast": 0.0 , #destroyed
                     } ,
                     "both": {
                         "CME": 3.618083294103339 , #destroyed
                         "slow": 1.5531765300203313e+02 , #destroyed
                         "fast": 2.901684756854498e+02 #destroyed
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 4.2163018750950135e+02 , #impact Sun
                         "slow": 4.2163018750950135e+02 , #impact Sun
                         "fast": 4.2163018750950135e+02 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #destroyed
                         "slow": 0.0 , #destroyed
                         "fast": 0.0 #destroyed
                     } ,
                     "both": {
                         "CME": 2.3508641524438683e+01 , #destroyed
                         "slow": 3.9897801034558944e+02 , #impact Sun
                         "fast": 4.12169761709957e+02 #impact Sun
                     }
                 }} }

#Slope calcs
# x1 = true_lifetime.get("particle2" , {}).get("size") 
# x2 = true_lifetime.get("particle3" , {}).get("size") 
# x3 = true_lifetime.get("particle1" , {}).get("size") 

# par2_sil_spu_slow = true_lifetime.get("particle2" , {}).get("silicate" , {}).get("sputtering" , {}).get("slow")
# par2_sil_spu_fast = true_lifetime.get("particle2" , {}).get("silicate" , {}).get("sputtering" , {}).get("fast")

# par2_car_spu_CME = true_lifetime.get("particle2" , {}).get("carbon" , {}).get("sputtering" , {}).get("CME")
# par2_car_spu_slow = true_lifetime.get("particle2" , {}).get("carbon" , {}).get("sputtering" , {}).get("slow")
# par2_car_spu_fast = true_lifetime.get("particle2" , {}).get("carbon" , {}).get("sputtering" , {}).get("fast")

# par3_sil_spu_slow = true_lifetime.get("particle3" , {}).get("silicate" , {}).get("sputtering" , {}).get("slow")
# par3_sil_spu_fast = true_lifetime.get("particle3" , {}).get("silicate" , {}).get("sputtering" , {}).get("fast")

# par3_car_spu_CME = true_lifetime.get("particle3" , {}).get("carbon" , {}).get("sputtering" , {}).get("CME")
# par3_car_spu_slow = true_lifetime.get("particle3" , {}).get("carbon" , {}).get("sputtering" , {}).get("slow")
# par3_car_spu_fast = true_lifetime.get("particle3" , {}).get("carbon" , {}).get("sputtering" , {}).get("fast")

# sil_spu_slow_slope , sil_spu_slow_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_sil_spu_slow , par3_sil_spu_slow)
# par1_sil_spu_slow = sil_spu_slow_slope * x3 + sil_spu_slow_intercept
    
# sil_spu_fast_slope , sil_spu_fast_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_sil_spu_fast , par3_sil_spu_fast)
# par1_sil_spu_fast = sil_spu_fast_slope * x3 + sil_spu_fast_intercept

# car_spu_CME_slope , car_spu_CME_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_car_spu_CME , par3_car_spu_CME)
# par1_car_spu_CME = car_spu_CME_slope * x3 + car_spu_CME_intercept

# car_spu_slow_slope , car_spu_slow_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_car_spu_slow , par3_car_spu_slow)
# par1_car_spu_slow = car_spu_slow_slope * x3 + car_spu_slow_intercept

# car_spu_fast_slope , car_spu_fast_intercept , _ , _ , _ = find_slope(x1 , x2 , par2_car_spu_fast , par3_car_spu_fast)
# par1_car_spu_fast = car_spu_fast_slope * x3 + car_spu_fast_intercept

if __name__ == "__main__":
    
    1
    

    