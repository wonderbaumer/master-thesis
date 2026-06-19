"""Numerical lifetimes in years, variable sputtering"""
pert_lifetime = {"C": {
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
                         "CME": 3226.745128539582 , #impact Sun
                         "slow": 3244.6933184534164 , #impact Sun
                         "fast": 3244.680954479293 , #impact Sun
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
                         "CME": 1514.511812171934 , #impact Sun
                         "slow": 1514.4249549140773 , #impact Sun
                         "fast": 1514.4249549140773 , #impact Sun
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
                         "CME": 11683.652027583568 , #impact Sun
                         "slow": 11682.719428376786 , #impact Sun
                         "fast": 11682.706624982955 , #impact Sun
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
                         "CME": 5919.0915016114695 , #impact Sun
                         "slow": 5919.0023714937915 , #impact Sun
                         "fast": 5919.001850265034 , #impact Sun
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
                         "CME": 202.310278880682 , #destroyed, singularity
                         "slow": 4053.3208468522434 , #destroyed, singularities
                         "fast": 4662.528128132865 , #destroyed, singularity
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
                         "CME": 3.008717951132692 , #destroyed
                         "slow": 121.56838112872384 , #destroyed
                         "fast": 243.2952861064858 , #destroyed
                     }
                 }}}