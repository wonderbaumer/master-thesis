import numpy as np

#Sputtering yields for astronomical silicate fast sw, slow sw, CME conditions, from Baumann, et al., 2020.
H_arr = np.array([0.0059 , 0.0034 , 0.0030 , 0.0209 , 0.0109 , 0.0047 , 0.0040 , 0.0346 , 0.0103 , 0.0049 , 0.0039 , 0.0305]) 
He_arr = np.array([0.0316 , 0.0179 , 0.0158 , 0.0997 , 0.0643 , 0.0354 , 0.0334 , 0.2199 , 0.0540 , 0.0281 , 0.0254 , 0.1696]) 
C_arr = np.array([0.1809 , 0.0999 , 0.0875 , 0.5699 , 0.3014 , 0.1603 , 0.1469 , 0.9641 , 0.2663 , 0.1440 , 0.1341 , 0.8738]) 
O_arr = np.array([0.2526 , 0.1401 , 0.1245 , 0.8138 , 0.4103 , 0.2171 , 0.1974 , 1.3400 , 0.3605 , 0.2008 , 0.1814 , 1.1800]) 
N_arr = np.array([0.2159 , 0.1183 , 0.1096 , 0.7044 , 0.3614 , 0.1924 , 0.1759 , 1.1600 , 0.3209 , 0.1748 , 0.1613 , 1.0400]) 
Fe_arr = np.array([0.8735 , 0.4798 , 0.4412 , 2.8600 , 0.9772 , 0.5234 , 0.4762 , 3.1800 , 0.9751 , 0.5217 , 0.4816 , 3.1700]) 
Ne_arr = np.array([0.3213 , 0.1784 , 0.1634 , 1.0600 , 0.4925 , 0.2704 , 0.2458 , 1.6000 , 0.4367 , 0.2343 , 0.2133 , 1.4100]) 
Mg_arr = np.array([0.3964 , 0.2142 , 0.1924 , 1.2800 , 0.5636 , 0.3010 , 0.2769 , 1.8200 , 0.5222 , 0.2875 , 0.2575 , 1.7100]) 
Si_arr = np.array([0.4783 , 0.2589 , 0.2393 , 1.5700 , 0.6416 , 0.3502 , 0.3141 , 2.0900 , 0.5772 , 0.3054 , 0.2834 , 1.8700]) 
S_arr = np.array([0.5407 , 0.2945 , 0.2688 , 1.7500 , 0.7136 , 0.3899 , 0.3431 , 2.3100 , 0.6676 , 0.3534 , 0.3217 , 2.1300]) 

#Elemental abundances in fast sw, slow sw and CME conditions from Killen, et al., 2012. All relative to hydrogen
n_tot = 1.0

n_He = np.array([0.02 , 0.04 , 0.3])
n_O = np.array([2.74e-4 , 4.86e-4 , 9.66e-3])
n_C = np.array([1.87e-4 , 3.19e-4 , 3.67e-3])
n_N = np.array([3.04e-5 , 4.19e-5 , 7.73e-4])
n_Ne = np.array([2.25e-5 , 4.95e-5 , 3.09e-3])
n_Mg = np.array([2.88e-5 , 6.81e-5 , 2.80e-3])
n_Si = np.array([3.15e-5 , 6.28e-5 , 1.74e-3])
n_S = np.array([1.53e-5 , 2.43e-5 , 1.18e-3])
n_Fe = np.array([2.52e-5 , 5.05e-5 , 7.05e-3])
n_H = n_tot - n_He - n_O - n_C - n_N - n_Ne - n_Mg - n_Si - n_S - n_Fe 

#Sputtering yields for carbon fast sw, slow sw and CME, from Baumann, et al., 2020.
H_C = np.array([0.0034 , 0.0076 , 0.0055])
He_C = np.array([0.0211 , 0.0567 , 0.0365])
C_C = np.array([0.1063 , 0.2389 , 0.1706])
O_C = np.array([0.1728 , 0.3614 , 0.2689])
N_C = np.array([0.1204 , 0.2606 , 0.1902])
Fe_C = np.array([1.3500 , 1.3500 , 1.3500])
Ne_C = np.array([0.3303 , 0.4620 , 0.4267])
Mg_C = np.array([0.6000 , 0.6000 , 0.6000])
Si_C = np.array([0.7200 , 0.7200 , 0.7200])
S_C = np.array([0.8300 , 0.8300 , 0.8300])

#Sputtering parameter, all solar wind conditions, all materials
sputter = {
    "silicate":{
        "fast": {
            "H": H_arr[0:4] * n_H[0] ,
            "He": He_arr[0:4] * n_He[0] ,
            "C": C_arr[0:4] * n_C[0] ,
            "O": O_arr[0:4] * n_O[0] ,
            "N": N_arr[0:4] * n_N[0] ,
            "Fe": Fe_arr[0:4] * n_Fe[0] ,
            "Ne": Ne_arr[0:4] * n_Ne[0] ,
            "Mg": Mg_arr[0:4] * n_Mg[0] ,
            "Si": Si_arr[0:4] * n_Si[0] ,
            "S": S_arr[0:4] * n_S[0] ,
        },
        "slow":{
            "H": H_arr[4:8] * n_H[1] ,
            "He": He_arr[4:8] * n_He[1] ,
            "C": C_arr[4:8] * n_C[1] ,
            "O": O_arr[4:8] * n_O[1] ,
            "N": N_arr[4:8] * n_N[1] ,
            "Fe": Fe_arr[4:8] * n_Fe[1] ,
            "Ne": Ne_arr[4:8] * n_Ne[1] ,
            "Mg": Mg_arr[4:8] * n_Mg[1] ,
            "Si": Si_arr[4:8] * n_Si[1] ,
            "S": S_arr[4:8] * n_S[1] ,
        },
        "CME":{
            "H": H_arr[8:12] * n_H[2] ,
            "He": He_arr[8:12] * n_He[2] ,
            "C": C_arr[8:12] * n_C[2] ,
            "O": O_arr[8:12] * n_O[2] ,
            "N": N_arr[8:12] * n_N[2] ,
            "Fe": Fe_arr[8:12] * n_Fe[2] ,
            "Ne": Ne_arr[8:12] * n_Ne[2] ,
            "Mg": Mg_arr[8:12] * n_Mg[2] ,
            "Si": Si_arr[8:12] * n_Si[2] ,
            "S": S_arr[8:12] * n_S[2] ,
        }},

    "carbon":{
        "fast":{
            "H": H_C[0] * n_H[0] ,
            "He": He_C[0] * n_He[0] ,
            "C": C_C[0] * n_C[0] ,
            "O": O_C[0] * n_O[0] ,
            "N": N_C[0] * n_N[0] ,
            "Fe": Fe_C[0] * n_Fe[0] ,
            "Ne": Ne_C[0] * n_Ne[0] ,
            "Mg": Mg_C[0] * n_Mg[0] ,
            "Si": Si_C[0] * n_Si[0] ,
            "S": S_C[0] * n_S[0] ,
            },

        "slow":{
            "H": H_C[1] * n_H[1] ,
            "He": He_C[1] * n_He[1] ,
            "C": C_C[1] * n_C[1] ,
            "O": O_C[1] * n_O[1] ,
            "N": N_C[1] * n_N[1] ,
            "Fe": Fe_C[1] * n_Fe[1] ,
            "Ne": Ne_C[1] * n_Ne[1] ,
            "Mg": Mg_C[1] * n_Mg[1] ,
            "Si": Si_C[1] * n_Si[1] ,
            "S": S_C[1] * n_S[1] ,
            },
        "CME":{
            "H": H_C[2] * n_H[2] ,
            "He": He_C[2] * n_He[2] ,
            "C": C_C[2] * n_C[2] ,
            "O": O_C[2] * n_O[2] ,
            "N": N_C[2] * n_N[2] ,
            "Fe": Fe_C[2] * n_Fe[2] ,
            "Ne": Ne_C[2] * n_Ne[2] ,
            "Mg": Mg_C[2] * n_Mg[2] ,
            "Si": Si_C[2] * n_Si[2] ,
            "S": S_C[2] * n_S[2] ,
            }
    }
}


