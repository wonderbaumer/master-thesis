import matplotlib.pyplot as plt
import numpy as np
from orbital_elements import eccentricity_sc, ecc_math
from config import t5 , t6 , t7 , c , V , B

def ecc_comps(x , y , vx , vy , beta , t , type):
    theta = np.atan2(y , x) #angle
    theta_cont = np.unwrap(theta) / (2 * np.pi) #continuous angle

    ecc_sc1 = eccentricity_sc(x , y , vx , vy , beta , t , type = type)

    ecc_m , time_arr , ecc_sc1 = ecc_math(x , y , t , ecc_sc1)
    ecc_m = ecc_m[::50]
    

    that = time_arr[::50]
    orbit = np.floor(theta_cont / (2*np.pi)).astype(int)
    theta_cont = theta_cont[::50]
    orbit = orbit[::50]

    ecc_sc1 = ecc_sc1[::50]
    
    plt.plot(theta_cont , ecc_sc1 * 10**5 , label = "Scaled e")
    plt.plot(theta_cont , ecc_m * 10**5, label = "Geometric e")
    
    plt.xlabel(r"$\hat{\theta}$ / $2\pi$")
    plt.ylabel(r"e value $\times 10^5$")
    plt.title("Scaled and geometric e as function of orbits")
    #plt.legend(bbox_to_anchor = (0.1 , 0.6))
    plt.legend()
    plt.show()

def test_plot(r_ssl , r_sfa , r_scme , r_csl , r_cfa , r_ccme , t):
    dt , t_tot = t
    time = np.arange(0 , t_tot , dt)

    plt.plot(time , r_ssl , label = "Slow solar wind" , linestyle = "-")
    plt.plot(time , r_sfa , label = "Fast solar wind" , linestyle = "-")
    plt.plot(time , r_scme , label = "CME" , linestyle = "-")

    plt.plot(time , r_csl , label = "Slow solar wind" , linestyle = "--")
    plt.plot(time , r_cfa , label = "Fast solar wind" , linestyle = "--")
    plt.plot(time , r_ccme , label = "CME" , linestyle = "--")

    plt.xlabel(r"$\hat{t}$")
    plt.ylabel(r"$\hat{r}$")

    plt.title("Numerical massloss and drag, silicate and carbon")
    plt.legend()
    plt.show()




if __name__ == "__main__":
    rk = np.load("Files/rk45_t5_masslossTrue_silicate_slow_noPR.npz")
    x1 , y1 , vx1 , vy1 , m1 , b1 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    
    s_f = np.load("Files/rk45_t5_masslossTrue_silicate_fast_noPR.npz")
    x2 , y2 , vx2 , vy2 , m2 , b2 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    
    s_cme = np.load("Files/rk45_t5_masslossTrue_silicate_CME_noPR.npz")
    x3 , y3 , vx3 , vy3 , m3 , b3 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]

    c_s = np.load("Files/rk45_t5_masslossTrue_carbon_slow_noPR.npz")
    x4 , y4 , vx4 , vy4 , m4 , b4 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    
    c_f = np.load("Files/rk45_t5_masslossTrue_carbon_fast_noPR.npz")
    x5 , y5 , vx5 , vy5 , m5 , b5 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    
    c_cme = np.load("Files/rk45_t5_masslossTrue_carbon_CME_noPR.npz")
    x6 , y6 , vx6 , vy6 , m6 , b6 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    
    r_ssl = np.sqrt(x1**2 + y1**2)
    r_sfa = np.sqrt(x2**2 + y2**2)
    r_scme = np.sqrt(x3**2 + y3**2)
    r_csl = np.sqrt(x4**2 + y4**2) 
    r_cfa = np.sqrt(x5**2 + y5**2)
    r_ccme = np.sqrt(x6**2 + y6**2)
    print(r_ssl)
    #test_plot(r_ssl , r_sfa , r_scme , r_csl , r_cfa , r_ccme , t5)
    