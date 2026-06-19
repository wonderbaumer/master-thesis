import numpy as np

"""Numerical lifetimes in years, constant sputtering"""
true_lifetime_csteps = {"A": {
                 "size": 1.54079 * 10**(-6) , 
                 "silicate": {
                     "pr": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 #placeholder
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 , #placeholder
                     } ,
                     "both": {
                         "CME": 3.22379014e+03 , #impact Sun
                         "slow": 3.21282706e+03 , #impact Sun
                         "fast": 3.21269320e+03 #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 #placeholder
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 #placeholder
                     } ,
                     "both": {
                         "CME": 1.50047145e+03 , #impact Sun
                         "slow": 1.49944971e+03 , #impact Sun
                         "fast": 1.49943621e+03 #impact Sun
                     }
                 }} , 
                 
                 "G": {
                     "size": 0.01220 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 #placeholder
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 #placeholder
                     } ,
                     "both": {
                         "CME": 206.6818304990219 , #impact Sun
                         "slow": 4806.665286067126 , #impact Sun
                         "fast": 4792.438426904735 #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 , #placeholder
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 , #placeholder
                     } ,
                     "both": {
                         "CME": 3.02165867 , #destroyed
                         "slow": 1.44143471e+02 , #impact Sun
                         "fast": 368.58088973612524 #impact Sun
                     }
                 }} , 
                 "H": {
                     "size": 0.00708 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 , #placeholder
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 , #placeholder
                     } ,
                     "both": {
                         "CME": 66.18015056856335 , #destroyed
                         "slow": 2840.9987086835345 , #destroyed
                         "fast": 5003.064814324503 #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 , #placeholder
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 , #placeholder
                     } ,
                     "both": {
                         "CME": 3.682629372272714 , #destroyed
                         "slow": 175.6578017656421 , #impact Sun
                         "fast": 421.6090358060843 #impact Sun
                     }
                 }} ,
                 "D": {
                     "size": 0.10165 * 10**(-6) ,
                 "silicate": {
                     "pr": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 , #placeholder
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 , #placeholder
                     } ,
                     "both": {
                         "CME": 73.17082948012936 , #impact Sun
                         "slow": 584.496817082005 , #impact Sun
                         "fast": 583.6234334919959 #impact Sun
                     }
                 } ,
                 "carbon": {
                     "pr": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 , #placeholder
                     } ,
                     "sputtering": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 , #placeholder
                     } ,
                     "both": {
                         "CME": 0.0 , #placeholder
                         "slow": 0.0 , #placeholder
                         "fast": 0.0 , #placeholder
                     }
                 }}}