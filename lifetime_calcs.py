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
                         "CME": 3.24501649e+03 , #impact Sun
                         "slow": 3.24501649e+03 , #impact Sun
                         "fast": 3.24501649e+03 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 2.36730114e+03 , #destroyed
                         "slow": 1.016238784334375e+05 , #assumed destroyed
                         "fast": 1.8985637134785356e+05 , #assumed destroyed
                     } ,
                     "both": {
                         "CME": 2.36730114e+03 , #destroyed
                         "slow": 3.20420150e+03 , #impact Sun
                         "fast": 3.22307505e+03 #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 1.51457811e+03 , #impact Sun
                         "slow": 1.51457811e+03 , #impact Sun
                         "fast": 1.51457811e+03 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 1.5381634133953323e+04 , #assumed destroyed
                         "slow": 7.335750199178321e+05 , #assumed destroyed
                         "fast": 1.8757798769500004e+06 #assumed destroyed
                     } ,
                     "both": {
                         "CME": 1.49295676e+03 , #impact Sun
                         "slow": 1.51405519e+03 , #impact Sun
                         "fast": 1.51437327e+03 #impact Sun
                     }
                 }} , 
                 
                 "particle2": {
                     "size": 0.01220 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 4.81684508e+03 , #impact Sun
                         "slow": 4.81684508e+03 , #impact Sun
                         "fast": 4.81684508e+03 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 1.72190836e+01 , #destroyed
                         "slow": 7.39183550e+02 , #destroyed
                         "fast": 1.38096192e+03 #destroyed
                     } ,
                     "both": {
                         "CME": 1.72190836e+01 , #destroyed
                         "slow": 7.39183550e+02 , #destroyed
                         "fast": 1.38096192e+03 #destroyed
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 4.10702687e+02 , #impact Sun
                         "slow": 4.10702687e+02 , #impact Sun
                         "fast": 4.10702687e+02 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 1.11881687e+02 , #destroyed
                         "slow": 5.33581866e+03 , #destroyed
                         "fast": 1.36438960e+04 #destroyed
                     } ,
                     "both": {
                         "CME": 1.05480156e+02 , #impact Sun
                         "slow": 3.41414628e+02 , #impact Sun
                         "fast": 3.77974424e+02 #impact Sun
                     }
                 }} , 
                 "particle3": {
                     "size": 0.00708 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 5.01579175e+03 , #impact Sun
                         "slow": 5.01579175e+03 , #impact Sun
                         "fast": 5.01579175e+03 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 9.34750253e+00 , #destroyed
                         "slow": 4.01271070e+02 , #destroyed
                         "fast": 7.49665041e+02 #destroyed
                     } ,
                     "both": {
                         "CME": 9.34750253e+00 , #destroyed
                         "slow": 4.01271070e+02 , #destroyed
                         "fast": 7.49665041e+02 #destroyed
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 4.21630188e+02 , #impact Sun
                         "slow": 4.21630188e+02 , #impact Sun
                         "fast": 4.21630188e+02 #impact Sun
                     } ,
                     "sputtering": {
                         "CME": 6.07357731e+01 , #destroyed
                         "slow": 2.89658727e+03 , #destroyed
                         "fast": 7.40668640e+03 #destroyed
                     } ,
                     "both": {
                         "CME": 6.07357731e+01 , #destroyed
                         "slow": 4.00647289e+02 , #impact Sun
                         "fast": 4.12456950e+02 #impact Sun
                     }
                 }}}

if __name__ == "__main__":
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

    print(car_spu_fast_intercept)