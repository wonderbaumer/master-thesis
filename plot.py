from matplotlib.pylab import beta
import matplotlib.pyplot as plt
import numpy as np
from constants import *
from config import *
import sys
sys.path.insert(1, 'C:/Users/Cecilie.Bamer/Documents/Project-paper/')
from analytical_functions import *
from energy import *
from forces import *

"""plotting params to make font sizes equal to font size in latex document"""
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

    plt.plot(b_init_vals[::-1] , epsilon[::-1]) #plots in reverse order
    plt.xlabel("B value")
    plt.ylabel("Epsilon value")
    plt.title("Epsilon vs B for size 500e-9 to 10e-6")
    plt.show()
    
"""comparing theta values between two solvers or one solver and perturbed expression"""
def ang_comps(solver1 , t , solver2 = None , theta_per = None):
    """input: solver1 (tuple), consisting of x and y coordinates for runge kutta 4 solver
              t (tuple), shape dt, t_tot for which x and y have been calculated
              solver2 (tuple), optional, default:None, x and y coordinates for Leapfrog solver
              theta_per (array), optional, default:None, theta values form perturbed expression
              
        returns: none   """
    
    x1 , y1 , vx1 , vy1 , m1 , b1 = solver1 #unpacking solver1
    theta1 = np.atan2(y1 , x1) #angle based on x and y
    theta1 = np.unwrap(theta1) #avoiding discontinuities in theta

    dt , t_tot = t #unpacking time params
    t = np.arange(0 , t_tot , dt) #time array

    """comparing theta as approximated by RK45 and Leapfrog solver"""
    if solver2 is not None:
        x2 , y2 , vx2 , vy2 , m2 , b2 = solver2 #solver2 unpacking
        theta2 = np.atan2(y2 , x2) #angle from x and y
        theta2 = np.unwrap(theta2) #removing discontinuities
        
        plt.plot(t / yr , theta1 , color = "blue", label = "RK45 solution")
        plt.plot(t / yr , theta2 , color = "red" , linestyle = "--" , label = "Leapfrog solution")
        plt.title("Theta, RK45 and Leapfrog comparison")

    """comparing theta calculated by RK45 to theta from perturbation expression"""
    if theta_per is not None:
        theta_per = np.unwrap(theta_per) #removing discontinuities
        plt.plot(t / yr , theta1 , color = "blue", label = "RK45 solution")
        plt.plot(t / yr, theta_per , color = "red" , linestyle = "--" , label = "Perturbed solution")
        plt.title("Theta, RK45 vs perturbed solution")

    plt.xlabel("Time (yr)")
    plt.ylabel("Angle (rad)")

    plt.legend()
    plt.show()

"""plotting radial distance as function of time, comparison of Leapfrog and RK45 solver, or RK45
and perturbation expression for r"""
def rad_comps(solver1 , t , solver2 = None , r_per = None):
    """input: solver1 (tuple), consisting of x and y, RK45 solver
              t (tuple), consisting of dt and t_tot, time over which x and y have been calculated
              solver2 (tuple), optional, default:None, x and y for Leapfrog solver
              r_per (array), optional, default:None, r from perturbed expression

        returns: none"""
    
    x1 , y1 , vx1 , vy1 , m1 , b1 = solver1 #solver1 unpacking
    r1 = np.sqrt(x1**2 + y1**2) #r calcs based on x and y

    dt , t_tot = t #t unpacking
    t = np.arange(0 , t_tot , dt) #time array 

    """comparing RK45 and Leapfrog approximations for r"""
    if solver2 is not None:
        x2 , y2 , vx2 , vy2 , m2 , b2 = solver2 #unpacking solver2
        r2 = np.sqrt(x2**2 + y2**2) #r calcs from x and y

        plt.plot(t / yr , r1 / au , color = "blue", label = "RK45 solution")
        plt.plot(t / yr , r2 / au , color = "red" , linestyle = "--" , label = "Leapfrog solution")
        plt.title("Radial distance, RK45 and Leapfrog comparisons")

    """comparing RK45 sol for r with r from perturbation expression"""
    if r_per is not None:
        plt.plot(t / yr , r1 / au , color = "blue" , label = "RK45 solution")
        plt.plot(t / yr , r_per , color = "red" , label = "Perturbed r")
        plt.title("Radial distance, numerical vs perturbed")

    plt.xlabel("Time (yr)")
    plt.ylabel("Radial distance (AU)")
    
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
    t = np.arange(0 , t_tot , dt) / yr #time array scaled to 1 year

    x , y , vx , vy , m , b = solver
    b = b / beta0 #RK45 beta hat
    b_per = b_per[: , 2] #first and second order beta hat perturbed expression
    b_analytical = b_analytical[: , 2] #first and second order beta hat analytical expression

    """comparing betahat from RK45 to betahat from perturbed and analytical expression"""
    if fw_err == False:
        plt.plot(t , b , color = "blue" , label = "RK45 beta")
        plt.plot(t , b_per , color = "red" , linestyle = "--" , label = "Perturbed beta")
        plt.plot(t , b_analytical , color = "orange" , linestyle = ":" , label = "Analytical beta")
        plt.xlabel("Time (yrs)")
        plt.ylabel("Beta hat")
        plt.title("Beta hat from RK45, perturbed and analytical expression")
    
    """comparing relative forward errors in beta hat from RK45-perturbed expression and 
    RK45-analytical expression"""
    if fw_err == True:
        fw_err_RK45_per = np.abs(solver - b_per) / np.abs(solver)
        fw_err_RK45_analytical = np.abs(solver - b_analytical) / np.abs(solver)
        plt.plot(t , fw_err_RK45_per , color = "blue" , label = "Rel fw error perturbed vs RK45")
        plt.plot(t , fw_err_RK45_analytical , color = "red" , linestyle = "--" , label = "Rel fw error analytical vs RK45")
        plt.xlabel("Time (yrs)")
        plt.ylabel("Relative forward error")
        plt.title("Relative forward error, perturbed vs RK45 sol and analytical vs RK45 sol")

    plt.legend()
    plt.show()

"""plots and compares energies between RK45 and Leapfrog solver"""
def energy_plot(t_arr , solver1 , solver2):
    """input: t_arr (tuple), consisting of dt and t_tot, time at which solver1 and 2 have been 
                             evaluated
              solver1 (tuple), consisting of x1, y1 , vx1 , vx2 , m1 , b_vals1, RK45 solver
              solver2 (tuple), x2 , y2 , vx1 , vx2 , m2 , b_vals2 , Leapfrog solver

       returns: none"""
    
    x1 , y1 , vx1 , vy1 , m1 , b_vals1 = solver1 #unpacking solver1
    x2 , y2 , vx2 , vy2 , m2 , b_vals2 = solver2 #unpacking solver2

    totenergy1 = tot_energy(x1 , y1 , vx1 , vy1 , m1 , b_vals1) #total energy calcs for solver1
    kinetic1 , potential1 = totenergy1 #unpacking into kinetic and potential energy
    tot1 = kinetic1 + potential1 #summing kinetic and potential energy into total energy

    totenergy2 = tot_energy(x2 , y2 , vx2 , vy2 , m2 , b_vals2) #total energy calcs solver2
    kinetic2 , potential2 = totenergy2 #unpacking into kinetic and potential energy
    tot2 = kinetic2 + potential2 #summing kinetic and potential energy into total energy
    
    dt , t_tot = t_arr #unpacking t_arr
    t_arr = np.arange(0 , t_tot , dt) / yr #time array scaled to 1 year
    
    #plots skipping 10 values for each iteration, for efficiency in plotting
    plt.plot(t_arr[::10] , kinetic1[::10] , label = "Kinetic RK45" , color = "blue" , linewidth = 2)
    plt.plot(t_arr[::10] , potential1[::10] , label = "Potential RK45" , color = "orange" , linewidth = 2)
    plt.plot(t_arr[::10] , tot1[::10] , label = "Total RK45" , color = "teal" , linewidth = 2)
    
    plt.plot(t_arr[::10] , kinetic2[::10] , label = "Kinetic Leapfrog" , color = "red" , linestyle = "--" , linewidth = 2)
    plt.plot(t_arr[::10] , potential2[::10] , label = "Potential Leapfrog" , color = "purple" , linestyle = "--" , linewidth = 2)
    plt.plot(t_arr[::10] , tot2[::10] , label = "Total Leapfrog" , color = "pink" , linestyle = "--" , linewidth = 2)
    
    plt.xlabel("Time (yr)")
    plt.ylabel("Energy (J)")
    plt.title("Runge Kutta 4(5) vs Leapfrog")
    plt.legend(loc = "upper right")

    plt.show()

if __name__ == "__main__":
    rk = np.load("C:/Users/Cecilie.Bamer/Documents/Project-paper/Files/rk45_500orbits_massloss.npz")
    x_r = rk["x"]
    y_r = rk["y"]
    vx_r = rk["vx"]
    vy_r = rk["vy"]
    m_r = rk["m"]
    b_r = rk["b"]

    lf = np.load("C:/Users/Cecilie.Bamer/Documents/Project-paper/Files/leapfrog_500orbits_massloss.npz")
    x_l = lf["x"]
    y_l = lf["y"]
    vx_l = lf["vx"]
    vy_l = lf["vy"]
    m_l = lf["m"]
    b_l = lf["b"]
    
    solver1 = b_r
    dt , t_tot = t4
    t_hat = np.arange(0 , t_tot , dt) / T
    b_per = betahat_pert(t_hat)
    b_analytical = betahat_analytical(t_hat)
    
    b_plot(solver1 , b_per , b_analytical , t4 , fw_err = True)

    
    
    
    
   
    
    
    
    


    
    

    