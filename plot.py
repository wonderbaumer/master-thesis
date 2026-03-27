import matplotlib.pyplot as plt
import numpy as np
from forces import beta
from config import init_cart , m_range , eps , t5 , t6 , t7 , T , B , r_betatest
from pert_functions import betahat_analytical , betahat_pert , rhat_pert , thetahat_pert , vrhat_pert ,omegahat_pert , perturbed_orbit , C0 , r , omega , theta , vr
from energy import tot_energy
from constants import dat_to_arr , sil_beta , car_beta , sputtering_lifetime , sputtering_yield , sw_flux , r_vals , M_ms , M_mc  
from forces_scaled import betahat , beta_real
from scipy.interpolate import PchipInterpolator as pchip

"""plotting params to adjust font sizes"""
plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
})

"""plots B values for masses corresponding to a range of initial particle sizes, 
and epsilon calculated same range of masses"""
def eps_init_beta():
    """input: none
       
       returns: none"""
    
    x , y = init_cart[0] , init_cart[1] #unpacking initial, unscaled cartesian values
    b_init_vals = beta(x , y , m_range)[1:] #calculating B values for mass range
    
    epsilon_slow = eps("silicate" , "slow" , "all" , m_range)[1:] #calculating epsilon for slow sw
    #epsilon_slow = epsilon_slow * 10**5 #scaling for better labelling

    epsilon_fast = eps("silicate" , "fast" , "all" , m_range)[1:] #calculating epsilon for fast sw
    #epsilon_fast = epsilon_fast * 10**5 #scaling for better labelling

    epsilon_cme = eps("silicate" , "CME" , "all" , m_range)[1:] #calculating epsilon for CME
    #epsilon_cme = epsilon_cme * 10**5 #scaling for better labelling

    epsilon_slowC = eps("carbon" , "slow" , "all" , m_range)[1:] #calculating epsilon for slow sw
    #epsilon_slowC = epsilon_slow * 10**5 #scaling for better labelling

    epsilon_fastC = eps("carbon" , "fast" , "all" , m_range)[1:] #calculating epsilon for fast sw
    #epsilon_fastC = epsilon_fast * 10**5 #scaling for better labelling

    epsilon_cmeC = eps("carbon" , "CME" , "all" , m_range)[1:] #calculating epsilon for CME
    #epsilon_cmeC = epsilon_cme * 10**5 #scaling for better labelling

    plt.plot(b_init_vals[::-1] , epsilon_slow[::-1] , color = "blue" , label = "Slow sw") #plots in reverse order
    plt.plot(b_init_vals[::-1] , epsilon_fast[::-1] , color = "red" , label = "Fast sw") #plots in reverse order
    plt.plot(b_init_vals[::-1] , epsilon_cme[::-1] , color = "green" , label = "CME") #plots in reverse order

    plt.plot(b_init_vals[::-1] , epsilon_slowC[::-1] , linestyle = "--" , color = "blue" ) #plots in reverse order
    plt.plot(b_init_vals[::-1] , epsilon_fastC[::-1] , linestyle = "--" , color = "red" ) #plots in reverse order
    plt.plot(b_init_vals[::-1] , epsilon_cmeC[::-1] , linestyle = "--" , color = "green") #plots in reverse order

    plt.xlabel(r"$B$")
    plt.ylabel(r"${\epsilon}$")
    plt.yscale("log")
    plt.title(r"${\epsilon}$ vs ${B}$, corresponding to size range $1~\mathrm{nm} \text{–} 50~\mu\mathrm{m}$, silicate and carbon")
    plt.legend(loc = "lower right")
    plt.show()
    
"""comparing thetahat values between RK4(5) and Leapfrog or RK4(5) and perturbed expression"""
def thetahat_comps(x1 , y1 , t , x2 = None , y2 = None , theta_per = None , species = None):
    """input: x1 (array), RK4(5) x vals
              y1 (array), RK4(5) y vals
              t (tuple), shape dt, t_tot for simulations
              x2 (array), optional, default:None, Leapfrog x vals
              y2 (array), optional, default:None, Leapfrog y vals
              theta_per (array), optional, default:None, thetahat values from perturbed expression
              
        returns: none   """

    theta1 = np.atan2(y1 , x1) #rk45 thetahat
    theta1 = np.unwrap(theta1) #avoiding discontinuities in thetahat
    
    dt , t_tot = t #unpacking time params
    t = np.arange(0 , t_tot , dt)  #t hat

    """comparing thetahat from RK4(5) and Leapfrog solver"""
    if x2 is not None and y2 is not None:
        theta2 = np.atan2(y2 , x2) #Leapfrog thetahat
        theta2 = np.unwrap(theta2) #removing discontinuities
        
        plt.plot(t[::10] , theta1[::10] , color = "blue", label = r"RK4(5) $\hat{\theta}$")
        plt.plot(t[::10] , theta2[::10] , color = "red" , linestyle = "--" , label = r"Leapfrog $\hat{\theta}$")
        plt.title(r"$\hat{\theta}$ from RK4(5) and Leapfrog solution")

    """comparing thetahat from RK4(5) perturbed thetahat"""
    if theta_per is not None:
        theta_per = np.unwrap(theta_per) #removing discontinuities
        plt.plot(t[::10] , theta1[::10] , color = "blue" , label = r"RK4(5) $\hat{\theta}$")
        plt.plot(t[::10] , theta_per[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\theta}$")
        plt.title(r"$\hat{\theta}$ from RK4(5) and perturbed solution")
    
    if species == "Silicate":
        plt.plot(t[::10] , theta1[::10])
        plt.title(r"$\hat{\theta}$ for silicate, real $\hat{\beta}$")
    
    if species == "Carbon":
        plt.plot(t[::10] , theta1[::10])
        plt.title(r"$\hat{\theta}$ for carbon, real $\hat{\beta}$")

    plt.xlabel("Number of orbits")
    plt.ylabel(r"$\hat{\theta}$")

    plt.legend()
    plt.show()

"""plotting rhat as function of orbits, comparison of Leapfrog and RK4(5) solver, or RK4(5)
and perturbed rhat"""
def rhat_comps(x1 , y1 , t , x2 = None , y2 = None , r_per = None , species = None):
    """input: x1 (array), RK4(5) x vals
              y1 (array), RK4(5) y vals
              t (tuple), consisting of dt and t_tot, time of simulations
              r_per (array), optional, perturbed rhat
              x2 (array), optional, default:None, Leapfrog x vals
              y2 (array), optional, default:None, Leapfrog y vals

        returns: none"""
    
    r1 = np.sqrt(x1**2 + y1**2) #r hat

    dt , t_tot = t #t unpacking
    t = np.arange(0 , t_tot , dt) #t hat

    orbit = round(len(t) / t_tot)
    orbit *=10

    """comparing RK4(5) and Leapfrog rhat"""
    if x2 is not None and y2 is not None:
        r2 = np.sqrt(x2**2 + y2**2) #r hat

        plt.plot(t[::10] , r1[::10] , color = "blue", label = r"RK4(5) $\hat{r}$")
        plt.plot(t[::10] , r2[::10] , color = "red" , linestyle = "--" , label = r"Leapfrog $\hat{r}$")
        plt.title(r"$\hat{r}$ from RK4(5) and Leapfrog solution")
        plt.legend()

    if r_per is not None: #comparing RK4(5) with perturbed rhat
        plt.plot(t[::10] , r1[::10] , color = "blue" , label = r"RK4(5) $\hat{r}$")
        plt.plot(t[::10] , r_per[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{r}$")
        plt.title(r"$\hat{r}$ from RK4(5) and perturbed solution")
        plt.legend()
    
    if species == "Silicate":
        r = np.sqrt(x1**2 + y1**2)
        plt.plot(t[::10] , r[::10])
        plt.title(r"$\hat{r}$ for silicate, real $\hat{\beta}$")
        
    
    if species == "Carbon":
        r = np.sqrt(x1**2 + y1**2)
        plt.plot(t[::10] , r[::10])
        plt.title(r"$\hat{r}$ for carbon, real $\hat{\beta}$")
        
    plt.xlabel(r"$\hat{t}$")
    plt.ylabel(r"$\hat{r}$")
    
    plt.show()

"plotting vhat from RK4(5) and perturbed expression, as function of t hat"
def vhat_comps(x , y , vx , vy , t , v_per):
    """input: x (array), RK4(5) x vals
              y (array), RK4(5) y vals
              vx (array), RK4(5) vx vals
              vy (array), RK4(5) vy vals
              t (tuple), consisting of dt , t_tot, time for simulations
              v_per (array), perturbed vhat

       returns: none"""
    
    theta_num = np.atan2(y , x) #thetahat
    theta_num = np.unwrap(theta_num) #avoiding discontinuities

    v_r = vx * np.cos(theta_num) + vy * np.sin(theta_num) #cartesian to radial vel
    #v_r = (x*vx + y*vy)/np.sqrt(x**2+y**2)
    vrplot = v_r * 10**5 #scaling for better labelling
    vperplot = v_per * 10**5 #scaling for better labelling
    #vrplot = np.sqrt(vx**2 + vy**2)
    dt , t_tot = t #time unpacking
    
    t = np.arange(0 , t_tot , dt)

    plt.plot(t[::10] , vrplot[::10] , color = "blue" , label = r"RK4(5) $\hat{v}$")
    plt.plot(t[::10] , vperplot[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{v}$")
    plt.xlabel("Number of orbits")
    plt.ylabel(r"$\hat{v} \times 10^5$")
    plt.title(r"RK4(5) and perturbed $\hat{v}$")
    plt.legend()
    plt.show()

"plotting omegahat from RK4(5) and perturbed expression as function of t hat"
def omegahat_comps(x , y , vx , vy , t , angvel):
    """input: x (array), RK4(5) x vals
              y (array), RK4(5) y vals
              vx (array), RK4(5) vx vals
              vy (array), RK4(5) vy vals
              t (tuple), consisting of dt, t_tot, time for simulations
              angvel (array), omegahat from perturbed expression

       returns: none"""
    
    theta_num = np.atan2(y , x) #thetahat RK4(5)
    theta_num = np.unwrap(theta_num) #avoiding discontinuities

    r = np.sqrt(x**2 + y**2) #r hat RK4(5)

    dt , t_tot = t #time unpacking

    t = np.arange(0 , t_tot , dt) #t hat
    angvel_num = (-vx * np.sin(theta_num) + vy * np.cos(theta_num)) / r    
    
    plt.plot(t[::10] , angvel_num[::10] , color = "blue" , label = "RK4(5)")
    plt.plot(t[::10] , angvel[::10] , color = "red" , linestyle = "--" , label = "Perturbed")
    plt.xlabel("Number of orbits")
    plt.ylabel(r"$\hat{\omega}$")
    plt.title(r"$\hat{\omega}$ RK4(5) vs perturbed solution")
    plt.legend()
    plt.show()

"""plotting betahat from RK4(5), perturbed and analytic expression.
Can compare betahat values or relative forward error RK4(5)-perturbed and RK4(5)-analytical"""
def b_plot(b_r , t , b_per , b_analytical , fw_err = False):
    """input: solver (.npz), RK4(5) solver file consisting of x, y, vx, vy, m, b
              b_per (array), betahat from perturbed expression
              b_analytical (array), betahat from analytical expression
              t (tuple), consisting of dt and t_tot at which solver params have been evaluated
              fw_err (bool), optional, default:False, user specifies if they want to plot relative 
                             forward errors RK4(5)-perturbed expression and RK4(5)-analytical expression
              
       returns: none"""
    
    dt , t_tot = t #unpacking t
    t = np.arange(0 , t_tot , dt) #that

    """comparing betahat from RK4(5) to betahat from perturbed and analytical expression"""
    if fw_err == False:
        plt.figure()
        plt.plot(t[::10] , b_r[::10] , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t[::10] , b_per[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ from RK4(5) and perturbed solution")
        plt.xlabel("Number of orbits")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.show()

        plt.figure()
        plt.plot(t[::10] , b_r[::10] , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t[::10] , b_analytical[::10] , color = "orange" , linestyle = "--" , label = r"Analytical $\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ from RK4(5) and analytical solution")
        plt.xlabel("Number of orbits")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.show()
        
        plt.figure()
        plt.plot(t[::10] , b_r[::10] , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t[::10] , b_per[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\beta}$")
        plt.plot(t[::10] , b_analytical[::10] , color = "orange" , linestyle = "--" , label = r"Analytical $\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ from RK4(5), perturbed and analytical solution")
        plt.xlabel("Number of orbits")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.show()

    
    """comparing relative forward errors in beta hat from RK45-perturbed expression and 
    RK45-analytical expression"""
    if fw_err == True:
        fw_err_RK45_per = np.abs(b_r - b_per) / np.abs(b_r)
        fw_err_RK45_analytical = np.abs(b_r - b_analytical) / np.abs(b_r)
        plt.plot(t[::10] , fw_err_RK45_per[::10] , color = "blue" , label = "Rel fw error RK4(5) vs perturbed")
        plt.plot(t[::10] , fw_err_RK45_analytical[::10] , color = "red" , linestyle = "--" , label = "Rel fw error RK4(5) vs analytical")
        plt.xlabel("Number of orbits")
        plt.ylabel(r"Relative forward error")
        plt.title("Relative forward error, RK4(5) vs perturbed and RK4(5) vs analytical")

        plt.legend()
        plt.show()

"""plots and compares energies between RK4(5) and Leapfrog solver"""
def energy_plot(solver1 , t , solver2 , fw_err = False):
    """input: solver1 (.npz), consisting of x1, y1 , vx1 , vx2 , m1 , b_vals1, RK4(5) solver
              t_arr (tuple), consisting of dt and t_tot, time at which solver1 and 2 have been 
                             evaluated
              solver2 (.npz), x2 , y2 , vx1 , vx2 , m2 , b_vals2 , Leapfrog solver

       returns: none"""
    
    x1 , y1 , vx1 , vy1 , m1 , b_vals1 = [solver1[k] for k in ("x","y","vx","vy","m","b")] #unpacking solver1
    x2 , y2 , vx2 , vy2 , m2 , b_vals2 = [solver2[k] for k in ("x","y","vx","vy","m","b")] #unpacking solver2
    x1 , y1 , vx1 , vy1 , m1 , b_vals1 = solver1
    x2 , y2 , vx2 , vy2 , m2 , b_vals2 = solver2
    totenergy1 = tot_energy(x1 , y1 , vx1 , vy1 , m1 , b_vals1) #total energy calcs for solver1
    kinetic1 , potential1 = totenergy1 #unpacking into kinetic and potential energy
    kinetic1 = kinetic1 
    potential1 = potential1  
    tot1 = kinetic1 + potential1 #summing kinetic and potential energy into total energy
    
    totenergy2 = tot_energy(x2 , y2 , vx2 , vy2 , m2 , b_vals2) #total energy calcs solver2
    kinetic2 , potential2 = totenergy2 #unpacking into kinetic and potential energy
    kinetic2 = kinetic2 
    potential2 = potential2 
    tot2 = kinetic2 + potential2 #summing kinetic and potential energy into total energy
    
    dt , t_tot = t #unpacking t_arr
    t = np.arange(0 , t_tot , dt) #that

    if fw_err == False:
        #plots skipping 10 values for each iteration, for efficiency in plotting
        plt.plot(t[::10] , kinetic1[::10] , label = "Kinetic RK4(5)" , color = "blue" , linewidth = 2)
        plt.plot(t[::10] , potential1[::10] , label = "Potential RK4(5)" , color = "orange" , linewidth = 2)
        plt.plot(t[::10] , tot1[::10] , label = "Total RK4(5)" , color = "teal" , linewidth = 2)
        plt.plot(t[::10] , kinetic2[::10] , label = "Kinetic Leapfrog" , color = "red" , linestyle = "--" , linewidth = 2)
        plt.plot(t[::10] , potential2[::10] , label = "Potential Leapfrog" , color = "purple" , linestyle = "--" , linewidth = 2)
        plt.plot(t[::10] , tot2[::10] , label = "Total Leapfrog" , color = "pink" , linestyle = "--" , linewidth = 2)
        plt.xlabel("Number of orbits")
        plt.ylabel("Energy")
        plt.title("RK4(5) vs Leapfrog")
        plt.legend(loc = "upper right" ,
               bbox_to_anchor = (1.0 , 0.8))
        plt.show()
    
    if fw_err == True: #plots RK4(5) sols and relative forward errors between RK4(5) and Leapfrog
        err_kin = np.abs(kinetic1 - kinetic2) / np.abs(kinetic1)
        err_pot = np.abs(potential1 - potential2) / np.abs(potential1)
        err_tot = np.abs(tot1 - tot2) / np.abs(tot1)

        plt.figure()
        plt.plot(t[::10] , kinetic1[::10] , label = "Kinetic RK4(5)" , color = "blue" , linewidth = 2)
        plt.plot(t[::10] , potential1[::10] , label = "Potential RK4(5)" , color = "orange" , linewidth = 2)
        plt.plot(t[::10] , tot1[::10] , label = "Total RK4(5)" , color = "teal" , linewidth = 2)
        plt.xlabel("Number of orbits")
        plt.ylabel("Energy")
        plt.title("RK4(5)")
        plt.legend(loc = "upper right" ,
               bbox_to_anchor = (1.0 , 0.9))
        plt.show()
        """
        plt.figure()
        plt.plot(t[::10] , err_kin[::10] , label = "Error kinetic" , color = "blue" , linewidth = 2)
        plt.plot(t[::10] , err_pot[::10]  , label = "Error potential" , color = "orange" , linewidth = 2)
        plt.plot(t[::10] , err_tot[::10] , label = "Error total energy" , color = "teal" , linewidth = 2)
        plt.xlabel("Number of orbits")
        plt.ylabel("Relative forward error")
        plt.title("Forward error, RK4(5) vs Leapfrog")
        plt.legend(loc = "upper right" ,
               bbox_to_anchor = (1.0 , 0.9))
        plt.show()
        """

"""plots beta curves for silicate and carbon"""
def beta_curves(interp = False , comp = False):
    """input: None
    
       returns: None"""
    
    sil_size , sil_betaval , _ = dat_to_arr(sil_beta) #fetching silicate size and beta values
    car_size , car_betaval , _ = dat_to_arr(car_beta) #fetching carbon size and beta values

    if interp == True: #compare interpolated function with true curve for silicate values
        interp = pchip(sil_size , sil_betaval)

        plt.plot(sil_size * 10**(-6) , sil_betaval , linestyle = "--" , label = "True curve")
        plt.plot(sil_size * 10**(-6) , interp(sil_size) , linestyle = "-" , label = "Interpolated function")
    
    if comp == True:

        plt.xscale("log")
        plt.yscale("log")
        plt.plot(sil_size * 10**(-6) , sil_betaval , color = "red" , linestyle = "-" , label = "Silicate")
        plt.plot(car_size * 10**(-6) , car_betaval , color = "blue" , linestyle = "--" , label = "Carbon")
    
    else:
        plt.xscale("log")
        plt.yscale("log")
        plt.scatter([1.54079 * 10**(-6) , 0.17508 * 10**(-6) , 0.04259 * 10**(-6)] , [0.1235 , 0.8560 , 0.2098] , c = "g")
        plt.scatter([1.54079 * 10**(-6) , 0.17508 * 10**(-6) , 0.04259 * 10**(-6)] , [0.2646 , 3.0589 , 1.6179] , c = "orange")
        plt.plot(sil_size * 10**(-6) , sil_betaval , color = "red" , linestyle = "-" , label = "Silicate")
        plt.plot(car_size * 10**(-6) , car_betaval , color = "blue" , linestyle = "--" , label = "Carbon")


    plt.title(r"$\beta$ versus particle size")
    plt.xlabel(r"Particle size (m)")
    plt.ylabel(r"$\beta$")
    plt.legend()
    plt.show()

def PR_spu_lifetime():
    """input: None
    
    returns: None"""
    
    f_sw , s_sw , CME_sw = sw_flux("fast") , sw_flux("slow") , sw_flux("CME")

    #Silicate
    sil_size , _ , sil_PR = dat_to_arr(sil_beta) #fetching silicate size and PR lifetime values
    fs_spu , ss_spu , CMEs_spu = sputtering_yield("silicate" , "fast") , sputtering_yield("silicate" , "slow") , sputtering_yield("silicate" , "CME")
    fs_lifetime , ss_lifetime , CMEs_lifetime = sputtering_lifetime(r_vals , f_sw , fs_spu , M_ms) , sputtering_lifetime(r_vals , s_sw , ss_spu , M_ms) , sputtering_lifetime(r_vals , CME_sw , CMEs_spu , M_ms)

    #Carbon
    car_size , _ , car_PR = dat_to_arr(car_beta) #fetching carbon size and PR lifetime values
    fc_spu , sc_spu , CMEc_spu = sputtering_yield("carbon" , "fast") , sputtering_yield("carbon" , "slow") , sputtering_yield("carbon" , "CME")
    fc_lifetime , sc_lifetime , CMEc_lifetime = sputtering_lifetime(r_vals , f_sw , fc_spu , M_mc) , sputtering_lifetime(r_vals , s_sw , sc_spu , M_mc) , sputtering_lifetime(r_vals , CME_sw , CMEc_spu , M_mc)

    plt.xscale("log")
    plt.yscale("log")
    plt.plot(sil_size * 10**(-6) , sil_PR , color = "black" , linestyle = "-" , label = "PR lifetime")
    plt.plot(car_size * 10**(-6) , car_PR , color = "black" , linestyle = "--")
    plt.plot(r_vals , fs_lifetime , color = "red" , linestyle = "-" , label = "Fast sw")
    plt.plot(r_vals , fc_lifetime , color = "red" , linestyle = "--")
    plt.plot(r_vals , ss_lifetime , color = "green" , linestyle = "-" , label = "Slow sw")
    plt.plot(r_vals , sc_lifetime , color = "green" , linestyle = "--")
    plt.plot(r_vals , CMEs_lifetime , color = "blue" , linestyle = "-" , label = "CME")
    plt.plot(r_vals , CMEc_lifetime , color = "blue" , linestyle = "--")

    plt.title(r"Poynting-Robertson and sputtering lifetimes")
    plt.xlabel(r"Particle size (m)")
    plt.ylabel(r"Lifetime (years)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    rk = np.load("Files/rk45_t6_silicateslowsw_realbeta.npz")
    x1 , y1 , vx1 , vy1 , m1 , b1 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    dt , t_tot = t6
    t = np.arange(0 , t_tot , dt)
    b_func = betahat_analytical(t)
    cst = C0(b_func)
    om = omega(t , b_func , cst)
    rad = r(t , b_func , cst , om)
    x = np.linspace(0 , 100 , 30)
    rhat_comps(x1 , y1 , t6 , species = "Silicate")
    
    
    


    

    
    
    
    
    