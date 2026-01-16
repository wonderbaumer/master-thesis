import matplotlib.pyplot as plt
import numpy as np
from constants import *
from config import *
import sys
sys.path.insert(1, 'C:/Users/cecil/Documents/Project-paper/')
from analytical_functions import *
from energy import *
from forces import *

"""plotting params to adjust font sizes"""
plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
})

"""plots initial beta values for a range of masses, and epsilon calculated for same range of masses"""
def eps_init_beta():
    """input: none
       
       returns: none"""
    
    x , y = init_cartesian[0] , init_cartesian[1] #unpacking initial cartesian values
    b_init_vals = beta(x , y , m_range) #calculating initial beta values for mass range

    epsilon = eps(m_range) #calculating epsilon
    epsilon = epsilon * 10**5 #scaling for better labelling

    plt.plot(b_init_vals[::-1] , epsilon[::-1]) #plots in reverse order
    plt.xlabel(r"$B$")
    plt.ylabel(r"$\epsilon \times 10^5$ ")
    plt.title(r"$\epsilon$ vs $B$, corresponding to size range $500^{-9} \text{–} 10^{-6}\,\mathrm{m}$")
    plt.show()
    
"""comparing theta values between two solvers or one solver and perturbed expression"""
def ang_comps(solver1 , t , solver2 = None , theta_per = None):
    """input: solver1 (tuple), consisting of x and y coordinates for runge kutta 4 solver
              t (tuple), shape dt, t_tot for which x and y have been calculated
              solver2 (tuple), optional, default:None, x and y coordinates for Leapfrog solver
              theta_per (array), optional, default:None, theta values form perturbed expression
              
        returns: none   """
    
    x1 , y1 , _ , _ , _ , _ = [solver1[k] for k in ("x","y","vx","vy","m","b")]  #unpacking solver1
    theta1 = np.atan2(y1 , x1) #angle based on x and y
    theta1 = np.unwrap(theta1) #avoiding discontinuities in theta
    
    dt , t_tot = t #unpacking time params
    t = np.arange(0 , t_tot , dt) / T  #t hat

    """comparing theta as approximated by RK45 and Leapfrog solver"""
    if solver2 is not None:
        x2 , y2 , _ , _ , _ , _ = [solver2[k] for k in ("x","y","vx","vy","m","b")]  #solver2 unpacking
        theta2 = np.atan2(y2 , x2) #angle from x and y
        theta2 = np.unwrap(theta2) #removing discontinuities
        
        plt.plot(t , theta1 , color = "blue", label = "RK4(5)")
        plt.plot(t , theta2 , color = "red" , linestyle = "--" , label = "Leapfrog")
        plt.title(r"$\hat{\theta}$, RK4(5) and Leapfrog")

    """comparing theta calculated by RK45 to theta from perturbation expression"""
    if theta_per is not None:
        theta_per = np.unwrap(theta_per) #removing discontinuities
        #plt.plot(t , theta1 , color = "blue", label = "RK4(5)")
        #plt.plot(t , theta_per, color = "red" , linestyle = "--" , label = "Perturbed solution")
        rel_fw_err = ((np.abs(theta1[1:]-theta_per[1:])) / np.abs(theta1[1:]))
        print(rel_fw_err)
        plt.plot(t[1:] , rel_fw_err)
        plt.title(r"Forward error, RK4(5) and perturbed $\hat{\theta}$")
        #plt.title(r"$\hat{\theta}$, RK4(5) vs perturbed solution")

    plt.xlabel("Number of orbits")
    #plt.ylabel(r"$\hat{\theta}$")
    #plt.yscale("log")
    plt.ylabel("Forward error")

    #plt.legend()
    plt.show()

"""plotting radial distance as function of time, comparison of Leapfrog and RK45 solver, or RK45
and perturbation expression for r"""
def rad_comps(solver1 , t , solver2 = None , r_per = None):
    """input: solver1 (tuple), consisting of x and y, RK45 solver
              t (tuple), consisting of dt and t_tot, time over which x and y have been calculated
              solver2 (tuple), optional, default:None, x and y for Leapfrog solver
              r_per (array), optional, default:None, r from perturbed expression

        returns: none"""
    
    x1 , y1 , _ , _ , _ , _ = [solver1[k] for k in ("x","y","vx","vy","m","b")] #solver1 unpacking
    r1 = np.sqrt(x1**2 + y1**2) #/ R #r hat
    dt , t_tot = t #t unpacking
    t = np.arange(0 , t_tot , dt) / T #t hat
    
    """comparing RK45 and Leapfrog approximations for r"""
    if solver2 is not None:
        x2 , y2 , _ , _ , _ , _ = [solver2[k] for k in ("x","y","vx","vy","m","b")] #unpacking solver2
        r2 = np.sqrt(x2**2 + y2**2) #/ R #r hat

        plt.plot(t , r1 , color = "blue", label = "RK4(5)")
        plt.plot(t , r2 , color = "red" , linestyle = "--" , label = "Leapfrog")
        plt.title(r"$\hat{r}$, RK4(5) and Leapfrog")

    else: #comparing RK45 sol for r with r from perturbation expression
        plt.plot(t , r1 , color = "blue" , label = "RK4(5)")
        plt.plot(t , r_per , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{r}$")
        #fw_err = (np.abs(r1-r_per) / np.abs(r1))
        #plt.plot(t , fw_err)
        #plt.title(r"$\hat{r}$, numerical vs perturbed")
        plt.title(r"RK4(5) and perturbed $\hat{r}$")
        
    plt.xlabel("Number of orbits")
    plt.ylabel(r"$\hat{r}$")
    #plt.ylabel("Forward error")
    
    plt.legend()
    plt.show()

"plotting radial velocity calculated by RK4(5) and perturbed expression, as function of t hat"
def rad_vel_comps(solver , t , v_per):
    """input: solver (.npz file), RK4(5) solutions
              t (tuple), consisting of dt,t_tot
              v_per (array), perturbed radial velocity

       returns: none"""
    
    x , y , vx , vy , _ , _ = [solver[k] for k in ("x","y","vx","vy","m","b")] #solver1 unpacking
    
    theta_num = np.atan2(y , x) #angle
    theta_num = np.unwrap(theta_num) #avoiding discontinuities

    v_r = vx * np.cos(theta_num) + vy * np.sin(theta_num) #cartesian to radial vel

    dt , t_tot = t #time unpacking

    t = np.arange(0 , t_tot , dt) / T #t hat

    plt.plot(t, v_r , color = "blue" , label = r"RK4(5) $\hat{v}$")
    plt.plot(t , v_per , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{v}$")
    plt.xlabel("Number of orbits")
    plt.ylabel(r"$\hat{v}$")
    plt.title(r"RK4(5) and perturbed $\hat{v}$")
    plt.legend()
    plt.show()

"plotting angular velocity calculated by RK4(5) and perturbed expression as function of t hat"
def ang_vel_comps(solver , t , angvel):
    """input: solver (.npz file), RK4(5) solutions
              t (tuple), consisting of dt, t_tot
              angvel (array), angular velocity perturbed expression

       returns: none"""
    
    x , y , vx , vy , _ , _ = [solver[k] for k in ("x","y","vx","vy","m","b")] #solver unpacking
    theta_num = np.atan2(y , x) #angle
    theta_num = np.unwrap(theta_num) #avoiding discontinuities
    v = np.sqrt(vx**2 + vy**2) #abs cartesian velocity
    r = np.sqrt(x**2 + y**2) #r hat

    dt , t_tot = t #time unpacking

    t = np.arange(0 , t_tot , dt) / T #t hat

    angvel_num = v / r #angular velocity
    
    plt.plot(t , angvel_num , label = "RK4(5)")
    plt.plot(t , angvel , label = "Perturbed")
    plt.xlabel("Number of orbits")
    plt.ylabel(r"$\hat{\omega}$")
    plt.title(r"$\hat{\omega}$ RK4(5) vs perturbed solution")
    plt.legend()
    plt.show()

"""plotting betahat as calculated numerically by RK45, from perturbed expression and analytic expression.
Can compare betahat values or relative forward error between RK45-perturbation and RK45-analytical"""
def b_plot(solver , b_per , b_analytical , t , fw_err = False):
    """input: solver (array), beta vals from RK45
              b_per (array), beta hat from perturbed expression
              b_analytical (array), beta hat from analytical expression
              t (tuple), consisting of dt and t_tot at which solver params have been evaluated
              fw_err (bool), optional, default:False, user specifies if they want to plot relative 
                             forward errors RK45-perturbed expression and RK45-analytical expression
              
       returns: none"""
    
    dt , t_tot = t #unpacking t
    t = np.arange(0 , t_tot , dt) / T #time array scaled to 1 orbital period
    _, _, _, _, _, b_r = [solver[k] for k in ("x","y","vx","vy","m","b")]

    #b = b_r / beta0 #RK45 beta hat
    b_per = b_per[: , 2] #first and second order beta hat perturbed expression

    """comparing betahat from RK4(5) to betahat from perturbed and analytical expression"""
    if fw_err == False:
        plt.figure()
        plt.plot(t , b_r , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t , b_per , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ from RK4(5) and perturbed expression")
        plt.xlabel("Number of orbits")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.show()

        plt.figure()
        plt.plot(t , b_r , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t , b_analytical , color = "orange" , linestyle = "--" , label = r"Analytical $\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ from RK4(5) and analytical expression")
        plt.xlabel("Number of orbits")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.show()
        
        plt.figure()
        plt.plot(t , b_r , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t , b_per , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\beta}$")
        plt.plot(t , b_analytical , color = "orange" , linestyle = "--" , label = r"Analytical $\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ from RK4(5), perturbed and analytical expression")
        plt.xlabel("Number of orbits")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.show()

    
    """comparing relative forward errors in beta hat from RK45-perturbed expression and 
    RK45-analytical expression"""
    if fw_err == True:
        fw_err_RK45_per = np.abs(b_r - b_per) / np.abs(b_r)
        fw_err_RK45_analytical = np.abs(b_r - b_analytical) / np.abs(b_r)
        plt.plot(t , fw_err_RK45_per , color = "blue" , label = "Rel fw error RK4(5) vs perturbed")
        plt.plot(t , fw_err_RK45_analytical , color = "red" , linestyle = "--" , label = "Rel fw error RK4(5) vs analytical")
        plt.xlabel("Number of orbits")
        plt.ylabel(r"Relative forward error")
        plt.title("Relative forward error, RK4(5) vs perturbed and RK4(5) vs analytical")

        plt.legend()
        plt.show()

"""plots and compares energies between RK45 and Leapfrog solver"""
def energy_plot(t_arr , solver1 , solver2 , fw_err = False):
    """input: t_arr (tuple), consisting of dt and t_tot, time at which solver1 and 2 have been 
                             evaluated
              solver1 (tuple), consisting of x1, y1 , vx1 , vx2 , m1 , b_vals1, RK45 solver
              solver2 (tuple), x2 , y2 , vx1 , vx2 , m2 , b_vals2 , Leapfrog solver

       returns: none"""
    
    x1 , y1 , vx1 , vy1 , m1 , b_vals1 = [solver1[k] for k in ("x","y","vx","vy","m","b")] #unpacking solver1
    x2 , y2 , vx2 , vy2 , m2 , b_vals2 = [solver2[k] for k in ("x","y","vx","vy","m","b")] #unpacking solver2

    totenergy1 = tot_energy(x1 , y1 , vx1 , vy1 , m1 , b_vals1) #total energy calcs for solver1
    kinetic1 , potential1 = totenergy1 #unpacking into kinetic and potential energy
    kinetic1 = kinetic1 #* 10**6
    potential1 = potential1 #* 10**6 
    tot1 = kinetic1 + potential1 #summing kinetic and potential energy into total energy

    totenergy2 = tot_energy(x2 , y2 , vx2 , vy2 , m2 , b_vals2) #total energy calcs solver2
    kinetic2 , potential2 = totenergy2 #unpacking into kinetic and potential energy
    kinetic2 = kinetic2 #* 10**6
    potential2 = potential2 #* 10**6
    tot2 = kinetic2 + potential2 #summing kinetic and potential energy into total energy
    dt , t_tot = t_arr #unpacking t_arr
    t_arr = np.arange(0 , t_tot , dt) / T #time array scaled to T

    print(tot1 , tot2)

    if fw_err == False:
        #plots skipping 10 values for each iteration, for efficiency in plotting
        plt.plot(t_arr[::10] , kinetic1[::10] , label = "Kinetic RK4(5)" , color = "blue" , linewidth = 2)
        plt.plot(t_arr[::10] , potential1[::10] , label = "Potential RK4(5)" , color = "orange" , linewidth = 2)
        plt.plot(t_arr[::10] , tot1[::10] , label = "Total RK4(5)" , color = "teal" , linewidth = 2)
    
        plt.plot(t_arr[::10] , kinetic2[::10] , label = "Kinetic Leapfrog" , color = "red" , linestyle = "--" , linewidth = 2)
        plt.plot(t_arr[::10] , potential2[::10] , label = "Potential Leapfrog" , color = "purple" , linestyle = "--" , linewidth = 2)
        plt.plot(t_arr[::10] , tot2[::10] , label = "Total Leapfrog" , color = "pink" , linestyle = "--" , linewidth = 2)
    
        plt.xlabel("Number of orbits")
        plt.ylabel("Energy")
        plt.title("RK4(5) vs Leapfrog")
        plt.legend(loc = "upper right" ,
               bbox_to_anchor = (1.0 , 0.8))
        #plt.show()
    
    if fw_err == True:
        err_kin = np.abs(kinetic1 - kinetic2) / np.abs(kinetic1)
        err_pot = np.abs(potential1 - potential1) / np.abs(potential1)
        err_tot = np.abs(tot1 - tot2) / np.abs(tot1)

        plt.figure()
        plt.plot(t_arr[::10] , kinetic1[::10] , label = "Kinetic RK4(5)" , color = "blue" , linewidth = 2)
        plt.plot(t_arr[::10] , potential1[::10] , label = "Potential RK4(5)" , color = "orange" , linewidth = 2)
        plt.plot(t_arr[::10] , tot1[::10] , label = "Total RK4(5)" , color = "teal" , linewidth = 2)
        plt.xlabel("Number of orbits")
        plt.ylabel("Energy")
        plt.title("RK4(5)")
        plt.legend(loc = "upper right" ,
               bbox_to_anchor = (1.0 , 0.9))
        plt.show()
        
        plt.figure()
        plt.plot(t_arr[::10] , err_kin[::10] , label = "Error kinetic" , color = "blue" , linewidth = 2)
        plt.plot(t_arr[::10] , err_pot[::10]  , label = "Error potential" , color = "orange" , linewidth = 2)
        plt.plot(t_arr[::10] , err_tot[::10] , label = "Error total energy" , color = "teal" , linewidth = 2)
        plt.xlabel("Number of orbits")
        plt.ylabel(r"Relative forward error $\times 10^7$")
        #plt.ylim(0 , 0.000001)
        plt.title("Forward error, RK4(5) vs Leapfrog")
        plt.legend(loc = "upper right" ,
               bbox_to_anchor = (1.0 , 0.9))
        plt.show()


if __name__ == "__main__":
    rk = np.load("C:/Users/cecil/Documents/Project-paper/Files/rk45_t7_masslossTrue_scaledeqs.npz")
    x_r = rk["x"]
    y_r = rk["y"]
    vx_r = rk["vx"]
    vy_r = rk["vy"]
    m_r = rk["m"]
    b_r = rk["b"]

    r_r = np.sqrt(x_r**2 + y_r**2) / R

    lf = np.load("C:/Users/cecil/Documents/Project-paper/Files/leapfrog_t6_masslossTrue_scaledeqs.npz")
    x_l = lf["x"]
    y_l = lf["y"]
    vx_l = lf["vx"]
    vy_l = lf["vy"]
    m_l = lf["m"]
    b_l = lf["b"]
    
    dt , t_tot = t7

    that = np.arange(0 , t_tot , dt) / T

    bper = betahat_pert(that)
    banalytical = betahat_analytical(that)

    theta_pert = angular_position(that)
    rper = radial_position(that)
    vper = rad_velocity(that)
    ang_vel = angular_vel(that)

    #b_plot(rk , bper , banalytical , t6 , fw_err = True)
    #rad_comps(rk , t6 , solver2 = None , r_per = rper)
    #energy_plot(t6 , rk , lf , fw_err = False)
    #ang_comps(rk , t7 , solver2 = None , theta_per = theta_pert)
    #rad_vel_comps(rk , t7 , vper)
    #ang_vel_comps(rk , t6 , ang_vel)



    

    
    
    
    
    