import matplotlib.pyplot as plt
import numpy as np
from forces import beta
from dust_properties import dust_properties
from pert_functions import perturbed_functions
from energy import tot_energy
from config import dat_to_arr , sil_beta , car_beta , r_vals , t5 , t6 , t7 , m_range , R , sil_size , sil_betaval , car_size , car_betaval , sil_PR , car_PR
from forces_scaled import betahat
from scipy.interpolate import PchipInterpolator as pchip
from polar_to_cart import polar_to_cartesian

"""plotting params to adjust font sizes"""
plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
})

"""plots B values for masses corresponding to a range of initial particle sizes, assuming beta=Frad/G, not real beta curve
and epsilon calculated same range of masses"""
def eps_init_beta(particle_obj):
    """input: none
       
       returns: none"""
    
    init_polar = np.array([R , 0 , 0 , particle_obj.V]) #initial polar coords
    init_cart = polar_to_cartesian(init_polar) #initial cart coords

    x , y = init_cart[0] , init_cart[1] #unpacking initial, unscaled cartesian values
    b_init_vals = beta(x , y , m_range)[1:] #calculating B values for mass range må fikse denne hvis real beta
    
    material = ["silicate" , "carbon"]
    sw_conds = ["slow" , "fast" , "CME"]
    eps_vals = []

    for m in material:
        for sw in sw_conds:
            eps = particle_obj.eps(m , sw)
            eps_vals.append({"material" : m ,
                             "sw_cond" : sw ,
                             "eps" : eps})
    
    for i in eps_vals:
        label = f"{i['material']} for {i['sw_cond']}"
        vals = i["eps"]
        plt.plot(b_init_vals[::10] , vals[::10] , label = label)

    plt.xlabel(r"$B$")
    plt.ylabel(r"${\epsilon}$")
    plt.yscale("log")
    plt.title(r"${\epsilon}$ vs ${B}$, corresponding to size range $1~\mathrm{nm} \text{–} 50~\mu\mathrm{m}$, silicate and carbon")
    plt.legend(loc = "lower right")
    plt.show()
    
"""comparing thetahat values between RK4(5) and Leapfrog or RK4(5) and perturbed expression"""
def thetahat_comps(x1 , y1 , t , x2 = None , y2 = None , particle_obj = None , material = None):
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
        plt.legend()

    """comparing thetahat from RK4(5) perturbed thetahat"""
    if particle_obj is not None:
        theta_per = particle_obj.theta()
        theta_per = np.unwrap(theta_per) #removing discontinuities
        plt.plot(t[::10] , theta1[::10] , color = "blue" , label = r"RK4(5) $\hat{\theta}$")
        plt.plot(t[::10] , theta_per[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\theta}$")
        plt.title(r"$\hat{\theta}$ from RK4(5) and perturbed solution")
        plt.legend()
    
    if material == "silicate":
        #t = t[:672160] #medium size, CME, impact Sun
        #t = t[:564630] #small size, CME, outside range
        #t = t[:2406880] #medium size, fast impact sun
        #t = t[:2395160] #t6 medium size, slow impact sun
        #t = t[:40535130] #t7 large, slow impact sun
        #t = t[:40868940] #t7 large, fast impact sun
        #t = t[:31398970] #t7 small, fast impact sun
        plt.plot(t[::10] , theta1[::10])
        plt.title(r"$\hat{\theta}$ for silicate, real $\hat{\beta}$")
    
    if material == "carbon":
        #t = t[:17637390] #t6 large fast impact sun
        #t = t[:17449500] #t6 large CME impact sun
        #t = t[:17634720] #t6 large slow impact sun  
        plt.plot(t[::10] , theta1[::10])
        plt.title(r"$\hat{\theta}$ for carbon, real $\hat{\beta}$")

    plt.xlabel(r"$\hat{t}$")
    plt.ylabel(r"$\hat{\theta}$")

    plt.show()

"""plotting rhat as function of orbits, comparison of Leapfrog and RK4(5) solver, or RK4(5)
and perturbed rhat"""
def rhat_comps(x1 , y1 , t , x2 = None , y2 = None , particle_obj = None , material = None):
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

    if particle_obj is not None: #comparing RK4(5) with perturbed rhat
        r_per , _ , _ = particle_obj.rad()
        rel_fw_err = np.abs(r1-r_per) / np.abs(r1)
        print(f"relative fw err:",rel_fw_err)
        plt.plot(t[::10] , r1[::10] , color = "blue" , label = r"RK4(5) $\hat{r}$")
        plt.plot(t[::10] , r_per[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{r}$")
        plt.title(r"$\hat{r}$ from RK4(5) and perturbed solution")
        plt.legend()
    
    if material == "silicate":
        r = np.sqrt(x1**2 + y1**2)
        #t = t[:672160] #medium size, CME, impact Sun
        #t = t[:564630] #small size, CME, outside range
        #t = t[:2406880] #medium size, fast impact sun
        #t = t[:2395160] #t6 medium size, slow impact sun
        #t = t[:40535130] #t7 large, slow impact sun
        #t = t[:40868940] #t7 large, fast impact sun
        #t = t[:31398970] #t7 small, fast impact sun
        plt.plot(t[::10] , r[::10])
        plt.title(r"$\hat{r}$ for silicate, real $\hat{\beta}$")
        
    
    if material == "carbon":
        #t = t[:17637390] #t6 large fast impact sun
        #t = t[:17449500] #t6 large CME impact sun
        #t = t[:17634720] #t6 large slow impact sun 
        r = np.sqrt(x1**2 + y1**2)
        plt.plot(t[::10] , r[::10])
        plt.title(r"$\hat{r}$ for carbon, real $\hat{\beta}$")
        
    plt.xlabel(r"$\hat{t}$")
    plt.ylabel(r"$\hat{r}$")
    
    plt.show()

"plotting vhat from RK4(5) and perturbed expression, as function of t hat"
def vhat_comps(x , y , t , vx , vy , particle_obj = None, material = None):
    """input: x (array), RK4(5) x vals
              y (array), RK4(5) y vals
              vx (array), RK4(5) vx vals
              vy (array), RK4(5) vy vals
              t (tuple), consisting of dt , t_tot, time for simulations
              v_per (array), perturbed vhat

       returns: none"""
    
    dt , t_tot = t #time unpacking
    
    t = np.arange(0 , t_tot , dt)

    theta_num = np.atan2(y , x) #thetahat
    theta_num = np.unwrap(theta_num) #avoiding discontinuities

    v_r = vx * np.cos(theta_num) + vy * np.sin(theta_num) #cartesian to radial vel
    
    if particle_obj is not None:
        v_per = particle_obj.vr()
    
        plt.plot(t[::10] , v_r[::10] , color = "blue" , label = r"RK4(5) $\hat{v}$")
        plt.plot(t[::10] , v_per[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{v}$")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v} \times 10^5$")
        plt.legend()
        plt.title(r"RK4(5) and perturbed $\hat{v}$")
    
    if material == "silicate":
        #t = t[:672160] #medium size, CME, impact Sun
        #t = t[:564630] #small size, CME, outside range
        #t = t[:2406880] #medium size, fast impact sun
        #t = t[:2395160] #t6 medium size, slow impact sun
        #t = t[:40535130] #t7 large, slow impact sun
        #t = t[:40868940] #t7 large, fast impact sun
        #t = t[:31398970] #t7 small, fast impact sun
        plt.plot(t[::10] , v_r[::10])
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_r$")
        plt.title(r"$\hat{v}_r$ for silicate, real $\hat{\beta}$")

    if material == "carbon":
        #t = t[:17637390] #t6 large fast impact sun
        #t = t[:17449500] #t6 large CME impact sun
        #t = t[:17634720] #t6 large slow impact sun 
        plt.plot(t[::10] , v_r[::10])
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_r$")
        plt.title(r"$\hat{v}_r$ for carbon, real $\hat{\beta}$")

    plt.show()

"plotting omegahat from RK4(5) and perturbed expression as function of t hat"
def omegahat_comps(x , y , vx , vy , t , particle_obj = None , material = None):
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

    if particle_obj is not None:
        angvel , _ , _ = particle_obj.omega()
        plt.plot(t[::10] , angvel_num[::10] , color = "blue" , label = "RK4(5)")
        plt.plot(t[::10] , angvel[::10] , color = "red" , linestyle = "--" , label = "Perturbed")
        plt.title(r"$\hat{\omega}$ RK4(5) vs perturbed solution")
        plt.legend()

    if material == "silicate":
        #t = t[:672160] #medium size, CME, impact Sun
        #t = t[:564630] #small size, CME, outside range
        #t = t[:2406880] #medium size, fast impact sun
        #t = t[:2395160] #t6 medium size, slow impact sun
        #t = t[:40535130] #t7 large, slow impact sun
        #t = t[:40868940] #t7 large, fast impact sun
        #t = t[:31398970] #t7 small, fast impact sun
        plt.plot(t[::10] , angvel_num[::10])
        plt.title(r"$\hat{\omega}$ for silicate, real $\hat{\beta}$")

    if material == "carbon":
        #t = t[:17637390] #t6 large fast impact sun
        #t = t[:17449500] #t6 large CME impact sun
        #t = t[:17634720] #t6 large slow impact sun 
        plt.plot(t[::10] , angvel_num[::10])
        plt.title(r"$\hat{\omega}$ for carbon, real $\hat{\beta}$")

    plt.xlabel(r"$\hat{t}$")
    plt.ylabel(r"$\hat{\omega}$")
    plt.show()

"""plotting betahat from RK4(5), perturbed and analytic expression.
Can compare betahat values or relative forward error RK4(5)-perturbed and RK4(5)-analytical"""
def b_plot(b_r , t , b_per = None , b_analytical = None, fw_err = False , material = None):
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
    if fw_err == False and b_per is not None and b_analytical is not None:
        plt.figure()
        plt.plot(t[::10] , b_r[::10] , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t[::10] , b_per[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ from RK4(5) and perturbed solution")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.show()

        plt.figure()
        plt.plot(t[::10] , b_r[::10] , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t[::10] , b_analytical[::10] , color = "orange" , linestyle = "--" , label = r"Analytical $\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ from RK4(5) and analytical solution")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.show()
        
        plt.figure()
        plt.plot(t[::10] , b_r[::10] , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t[::10] , b_per[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\beta}$")
        plt.plot(t[::10] , b_analytical[::10] , color = "orange" , linestyle = "--" , label = r"Analytical $\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ from RK4(5), perturbed and analytical solution")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.show()
    
    """comparing relative forward errors in beta hat from RK45-perturbed expression and 
    RK45-analytical expression"""
    if fw_err == True and b_per is not None and b_analytical is not None:
        fw_err_RK45_per = np.abs(b_r - b_per) / np.abs(b_r)
        fw_err_RK45_analytical = np.abs(b_r - b_analytical) / np.abs(b_r)
        plt.plot(t[::10] , fw_err_RK45_per[::10] , color = "blue" , label = "Rel fw error RK4(5) vs perturbed")
        plt.plot(t[::10] , fw_err_RK45_analytical[::10] , color = "red" , linestyle = "--" , label = "Rel fw error RK4(5) vs analytical")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"Relative forward error")
        plt.title("Relative forward error, RK4(5) vs perturbed and RK4(5) vs analytical")

        plt.legend()
        plt.show()
    
    if material == "silicate":
        #t = t[:672160] #t5 medium size, CME, impact Sun
        #t = t[:564630] #t5 small size, CME outside range
        #t = t[:2406880] #t6 medium size, fast impact sun
        #t = t[:2395160] #t6 medium size, slow impact sun
        #t = t[:40535130] #t7 large, slow impact sun
        #t = t[:40868940] #t7 large, fast impact sun
        #t = t[:31398970] #t7 small, fast impact sun
        plt.plot(t[::10] , b_r[::10])
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ silicate")
        plt.show()
    
    if material == "carbon":
        #t = t[:17637390] #t6 large fast impact sun
        #t = t[:17449500] #t6 large CME impact sun
        #t = t[:17634720] #t6 large slow impact sun 
        plt.plot(t[::10] , b_r[::10])
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ carbon")
        plt.show()

"""plots and compares energies between RK4(5) and Leapfrog solver"""
def energy_plot(solver1 , t , solver2 , particle_obj , fw_err = False):
    """input: solver1 (.npz), consisting of x1, y1 , vx1 , vx2 , m1 , b_vals1, RK4(5) solver
              t_arr (tuple), consisting of dt and t_tot, time at which solver1 and 2 have been 
                             evaluated
              solver2 (.npz), x2 , y2 , vx1 , vx2 , m2 , b_vals2 , Leapfrog solver

       returns: none"""
    
    x1 , y1 , vx1 , vy1 , m1 , b_vals1 = [solver1[k] for k in ("x","y","vx","vy","m","b")] #unpacking solver1
    x2 , y2 , vx2 , vy2 , m2 , b_vals2 = [solver2[k] for k in ("x","y","vx","vy","m","b")] #unpacking solver2
    x1 , y1 , vx1 , vy1 , m1 , b_vals1 = solver1
    x2 , y2 , vx2 , vy2 , m2 , b_vals2 = solver2
    totenergy1 = tot_energy(x1 , y1 , vx1 , vy1 , m1 , b_vals1 , particle_obj) #total energy calcs for solver1
    kinetic1 , potential1 = totenergy1 #unpacking into kinetic and potential energy
    kinetic1 = kinetic1 
    potential1 = potential1  
    tot1 = kinetic1 + potential1 #summing kinetic and potential energy into total energy
    
    totenergy2 = tot_energy(x2 , y2 , vx2 , vy2 , m2 , b_vals2 , particle_obj) #total energy calcs solver2
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
        plt.plot(t[::10] , err_kin[::10] , label = "Error kinetic" , color = "blue" , linewidth = 2)
        plt.plot(t[::10] , err_pot[::10]  , label = "Error potential" , color = "orange" , linewidth = 2)
        plt.plot(t[::10] , err_tot[::10] , label = "Error total energy" , color = "teal" , linewidth = 2)
        plt.xlabel("Number of orbits")
        plt.ylabel("Relative forward error")
        plt.title("Forward error, RK4(5) vs Leapfrog")
        plt.legend(loc = "upper right" ,
               bbox_to_anchor = (1.0 , 0.9))
        plt.show()
        
"""plots beta curves for silicate and carbon"""
def beta_curves(interp = False , comp = False , sample_pts = False):
    """input: None
    
       returns: None"""
    

    """compare interpolated function with true curve for silicate values"""
    if interp == True: 
        interp = pchip(sil_size , sil_betaval)

        plt.plot(sil_size * 10**(-6) , sil_betaval , linestyle = "--" , label = "True curve")
        plt.plot(sil_size * 10**(-6) , interp(sil_size) , linestyle = "-" , label = "Interpolated function")
    
    """comparing real silicate and carbon beta curve"""
    if comp == True:

        plt.xscale("log")
        plt.yscale("log")
        plt.plot(sil_size * 10**(-6) , sil_betaval , color = "red" , linestyle = "-" , label = "Silicate")
        plt.plot(car_size * 10**(-6) , car_betaval , color = "blue" , linestyle = "--" , label = "Carbon")
    
    """curve with sampling points"""
    if sample_pts == True:
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

def PR_spu_lifetime(particle_obj):
    """input: None
    
    returns: None"""

    material = ["silicate" , "carbon"]
    sw_conds = ["slow" , "fast" , "CME"]
    tsp_vals = []
    ls = {"silicate" : "-" , 
          "carbon" : "--"}
    cl = {"slow" : "green" ,
          "fast" : "red" ,
          "CME" : "blue"}

    for m in material:
        for sw in sw_conds:
            t_sp = particle_obj.sputtering_lifetime(m , sw)
            tsp_vals.append({"material" : m ,
                             "sw_cond" : sw ,
                             "tsp" : t_sp})
    for i in tsp_vals:
        label = f"{i['material']} for {i['sw_cond']}"
        vals = i["tsp"]
        plt.plot(r_vals[::10] , vals[::10] , color = cl[i['sw']] , linestyle = ls[i['material']] , label = label)


    plt.xscale("log")
    plt.yscale("log")
    plt.plot(sil_size * 10**(-6) , sil_PR , color = "black" , linestyle = "-" , label = "PR lifetime")
    plt.plot(car_size * 10**(-6) , car_PR , color = "black" , linestyle = "--")


    plt.title(r"Poynting-Robertson and sputtering lifetimes")
    plt.xlabel(r"Particle size (m)")
    plt.ylabel(r"Lifetime (years)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    res = np.load("Files/rk45_t7_small_silicate_slowsw.npz")
    x , y , vx , vy , m , b = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
    #thetahat_comps(x , y , t7 , material = "silicate")
    #omegahat_comps(x , y , vx , vy , t7 , material = "silicate")
    #vhat_comps(x , y , t7 , vx , vy , material = "silicate")
    #rhat_comps(x , y , t7 , material = "silicate")
    b_plot(b , t7 , material = "silicate")
    
    # par = dust_properties("silicate" , "slow" , "all" , "large")
    # dt , t_tot = t5
    # t = np.arange(0 , t_tot , dt)
    
    # pert = perturbed_functions(drag_coeff = 13.238542083399384 , particle = par , t = t5)
    # #omegahat_comps(x , y , vx , vy , t5 , pert)
    # rhat_comps(x , y , t5 , particle_obj = pert)
    
    
    
    


    

    
    
    
    
    