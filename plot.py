import matplotlib.pyplot as plt
import numpy as np
from forces import beta
from dust_properties import dust_properties
from pert_variable_eps import perturbed_functions
from energy import tot_energy
from config import car_betaval_bound , car_size_bound , init_vals , sil_beta , car_beta , t5 , t6 , t7 , sil_size , sil_betaval , car_size , car_betaval , sil_PR , car_PR , sil_mass , car_mass , tau_car , tau_sil
from forces_scaled import betahat
from scipy.interpolate import PchipInterpolator as pchip
from polar_to_cart import polar_to_cartesian
import os
from scipy.constants import G
from lifetime_calcs import true_lifetime , true_lifetime_variableeps
import matplotlib.patches as mpatches
from orbital_elements import ecc_calcs , ecc_scaled

"""plotting params to adjust font sizes"""
plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
})

"""plots B values for masses corresponding to a range of initial particle sizes, assuming real beta curve
and epsilon calculated same range of masses"""
def eps_init_betareal():
    """input: none
       
       returns: none"""
    
    combs = {"silicate" : {
        "fast" : {"linestyle" : "-" , "color" : "red"} ,
        "slow" : {"linestyle" : "-" , "color" : "green"} ,
        "CME" :  {"linestyle" : "-" , "color" : "blue"}} ,

            "carbon" : {
        "fast" : {"linestyle" : "--" , "color" : "red"} ,
        "slow" : {"linestyle" : "--" , "color" : "green"} ,
        "CME" :  {"linestyle" : "--" , "color" : "blue"}}}
    
    
    for mat , sw in combs.items():

        for sw_cond , styles in sw.items():
            size_ranges = (sil_size , sil_betaval) if mat == "silicate" else (car_size_bound , car_betaval_bound)
            par = dust_properties(mat , sw_cond , size = None , size_range = size_ranges)
            init_beta = par.B
            epsilon = par.eps()
            delta = par.delta
            
            
            plt.plot(init_beta , epsilon , color = styles["color"] , linestyle = styles["linestyle"])
            plt.plot(init_beta , par.delta , color = "purple" , linestyle = styles["linestyle"])
            

        purple_patch = mpatches.Patch(color = "purple" , label = r"${\delta}$")
        blue_patch = mpatches.Patch(color = "blue" , label = r"${\epsilon_0}$ CME")
        red_patch = mpatches.Patch(color = "red" , label = r"${\epsilon_0}$ Fast")
        green_patch = mpatches.Patch(color = "green" , label = r"${\epsilon_0}$ Slow")

        handles, labels = plt.gca().get_legend_handles_labels()
        handles.extend([purple_patch , blue_patch , red_patch , green_patch])    
        plt.loglog()        
        plt.title(fr"${{\epsilon_0}}$ and ${{\delta}}$ vs B for {mat}")
        plt.xlabel("B")
        plt.ylabel("Value")
        plt.ylim(10**(-9) , 1)
        plt.legend(handles = handles)
        plt.savefig(f"Plots/{mat}_epsilonvsbeta.png", dpi = 300 , bbox_inches = 'tight')
        plt.show()
    
"""comparing thetahat values between RK4(5) and Leapfrog or RK4(5) and perturbed expression"""
def thetahat_comps(file_path , file_path_comp = None , pert = None , material = None):
    """input: x1 (array), RK4(5) x vals
              y1 (array), RK4(5) y vals
              t (tuple), shape dt, t_tot for simulations
              x2 (array), optional, default:None, Leapfrog x vals
              y2 (array), optional, default:None, Leapfrog y vals
              theta_per (array), optional, default:None, thetahat values from perturbed expression
              
        returns: none   """
    
    res = np.load(file_path)
    x , y , _ , _ , _ , _ , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b" , "t")]
    r = np.sqrt(x**2 + y**2)

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_theta.png"

    theta1 = np.atan2(y , x) #rk45 thetahat
    theta1 = np.unwrap(theta1) #avoiding discontinuities in thetahat
    
    """comparing thetahat from RK4(5) and Leapfrog solver"""
    if file_path_comp is not None:
        base_name1 = os.path.splitext(os.path.basename(file_path_comp))[0]
        save_path1 = f"Plots/{base_name}_theta_vs_{base_name1}.png"
        res1 = np.load(file_path_comp)
        x1 , y1 , _ , _ , _ , _ = [res1[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
        theta2 = np.atan2(y1 , x1) #Leapfrog thetahat
        theta2 = np.unwrap(theta2) #removing discontinuities
        
        plt.plot(t[r >= 0.1] , theta1[r >= 0.1] , color = "blue", label = r"RK4(5) $\hat{\theta}$")
        plt.plot(t[r >= 0.1] , theta2[r >= 0.1] , color = "red" , linestyle = "--" , label = r"Leapfrog $\hat{\theta}$")
        plt.title(r"$\hat{\theta}$ from RK4(5) and Leapfrog solution")
        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')

    """comparing thetahat from RK4(5) perturbed thetahat"""
    if pert is not None:
        theta_per = pert
        save_path2 = f"Plots/{base_name}_theta_vs_perturbed.png"
        theta_per = np.unwrap(theta_per) #removing discontinuities

        plt.plot(t[r >= 0.1] , theta1[r >= 0.1] , color = "blue" , label = r"RK4(5) $\hat{\theta}$")
        plt.plot(t[r >= 0.1] , theta_per[r >= 0.1] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\theta}$")
        
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\theta}$")
        plt.title(rf"$\hat{{\theta}}$ from RK4(5) and perturbed solution, {material}")
        plt.legend()
        # plt.savefig(save_path2 , dpi = 300 , bbox_inches = 'tight')
    
    if file_path_comp is None and pert is None:
        # plt.plot(t[r >= 0.1] , theta1[r >= 0.1] , color = "blue")
        plt.plot(t , theta1 , color = "blue")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\theta}$")

        plt.title(rf"$\hat{{\theta}}$ for {material}")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    plt.show()

"""plotting rhat as function of orbits, comparison of Leapfrog and RK4(5) solver, or RK4(5)
and perturbed rhat"""
def rhat_comps(file_path , material , file_path_comp = None , pert = None):
    """input: x1 (array), RK4(5) x vals
              y1 (array), RK4(5) y vals
              t (tuple), consisting of dt and t_tot, time of simulations
              r_per (array), optional, perturbed rhat
              x2 (array), optional, default:None, Leapfrog x vals
              y2 (array), optional, default:None, Leapfrog y vals

        returns: none"""
    
    res = np.load(file_path)
    x , y , _ , _ , _ , _ , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    r = np.sqrt(x**2 + y**2) #r hat
    
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_r.png"

    if pert is not None and file_path_comp is None: #comparing RK4(5) with perturbed rhat
        r_per = pert

        rel_fw_err = np.abs(r - r_per) / np.abs(r)
        # print(rel_fw_err)
        save_path2 = f"Plots/{base_name}_r_vs_perturbed.png"
        
        plt.plot(t[r >= 0.1] , r[r >= 0.1] , color = "blue" , label = r"RK4(5) $\hat{r}$")
        plt.plot(t[r >= 0.1] , r_per[r >= 0.1] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{r}$")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{r}$")
        plt.title(fr"{material.capitalize()} $\hat{{r}}$, RK4(5) and perturbed solution")
        plt.legend()
        plt.savefig(save_path2 , dpi = 300 , bbox_inches = 'tight')
    
    if file_path_comp is None and pert is None:
        # plt.plot(t[r >= 0.1] , r[r >= 0.1] , color = "blue" , label = fr"$\hat{{r}}$ for {material}")
        plt.plot(t , r , color = "blue" , label = fr"$\hat{{r}}$ for {material}")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{r}$")
        plt.title(fr"$\hat{{r}}$ for {material}")
        
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
    
    plt.show()

"plotting vhat from RK4(5) and perturbed expression, as function of t hat"
def vhat_comps(file_path , pert = None , material = None):
    """input: x (array), RK4(5) x vals
              y (array), RK4(5) y vals
              vx (array), RK4(5) vx vals
              vy (array), RK4(5) vy vals
              t (tuple), consisting of dt , t_tot, time for simulations
              v_per (array), perturbed vhat

       returns: none"""

    res = np.load(file_path)
    x , y , vx , vy , _ , _ , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    r = np.sqrt(x**2 + y**2)

    theta_num = np.atan2(y , x) #thetahat
    theta_num = np.unwrap(theta_num) #avoiding discontinuities

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_vr.png"

    v_r = vx * np.cos(theta_num) + vy * np.sin(theta_num) #cartesian to radial vel
    
    if pert is not None:
        v_per = pert
        # rel_fw_err = np.abs(v_r - v_per) / np.abs(v_r)
        save_path1 = f"Plots/{base_name}_vr_vs_perturbed.png"
        plt.plot(t[r >= 0.1] , v_r[r >= 0.1] , color = "blue" , label = r"RK4(5) $\hat{v}_r$")
        plt.plot(t[r >= 0.1] , v_per[r >= 0.1] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{v}_r$")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_r$")

        plt.legend()
        plt.title(rf"{material.capitalize()} $\hat{{v}}_r$, RK4(5) and perturbed solution")
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')
    
    else:
        # plt.plot(t[r >= 0.1] , v_r[r >= 0.1] , color = "blue")
        plt.plot(t , v_r , color = "blue")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_r$")
        plt.title(fr"$\hat{{v}}_r$ for {material}")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    # plt.show()

"plotting omegahat from RK4(5) and perturbed expression as function of t hat"
def omegahat_comps(file_path , pert = None , material = None):
    """input: x (array), RK4(5) x vals
              y (array), RK4(5) y vals
              vx (array), RK4(5) vx vals
              vy (array), RK4(5) vy vals
              t (tuple), consisting of dt, t_tot, time for simulations
              angvel (array), omegahat from perturbed expression

       returns: none"""
    
    res = np.load(file_path)
    x , y , vx , vy , _ , _ , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    r = np.sqrt(x**2 + y**2)

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_omega.png"
    
    theta_num = np.atan2(y , x) #thetahat RK4(5)
    theta_num = np.unwrap(theta_num) #avoiding discontinuities

    angvel_num = (-vx * np.sin(theta_num) + vy * np.cos(theta_num)) / r    

    if pert is not None:
        angvel = pert
        plt.plot(t[r >= 0.1] , angvel_num[r >= 0.1] , color = "blue" , label = r"RK4(5) $\hat{\omega}$")
        plt.plot(t[r >= 0.1] , angvel[r >= 0.1] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\omega}$")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\omega}$")
        plt.title(rf"{material.capitalize()} $\hat{{\omega}}$, RK4(5) and perturbed solution")

        plt.legend()
        save_path = f"Plots/{base_name}_omega_vs_perturbed.png"
        
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    if pert is None:
        # plt.plot(t[r >= 0.1] , angvel_num[r >= 0.1] , color = "blue")
        plt.plot(t , angvel_num , color = "blue")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\omega}$")
        
        plt.title(rf"$\hat{{\omega}}$ for {material}")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    plt.show()

"""plotting betahat from RK4(5), perturbed and analytic expression.
Can compare betahat values or relative forward error RK4(5)-perturbed and RK4(5)-analytical"""
def b_plot(file_path , b_per = None , b_analytical = None, fw_err = False , material = None):
    """input: solver (.npz), RK4(5) solver file consisting of x, y, vx, vy, m, b
              b_per (array), betahat from perturbed expression
              b_analytical (array), betahat from analytical expression
              t (tuple), consisting of dt and t_tot at which solver params have been evaluated
              fw_err (bool), optional, default:False, user specifies if they want to plot relative 
                             forward errors RK4(5)-perturbed expression and RK4(5)-analytical expression
              
       returns: none"""
    
    res = np.load(file_path)
    x , y , _ , _ , _ , b , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    r = np.sqrt(x**2 + y**2)

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path2 = f"Plots/{base_name}_beta.png"

    if fw_err == False and b_per is not None and b_analytical is None:
        b_pert = b_per

        save_path3 = f"Plots/{base_name}_beta_interp_vs_expression.png"
        plt.plot(t[r >= 0.1] , b_pert[r >= 0.1] , color = "red" , label = r"$\hat{\beta}$ expression")
        plt.plot(t[r >= 0.1] , b[r >= 0.1] , color = "blue" , linestyle = "--" , label = r"Interpolated $\hat{\beta}$ curve")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.title(fr"$\hat{{\beta}}$ from interpolated function and from expression, {material}")
        plt.legend()
        plt.savefig(save_path3 , dpi = 300 , bbox_inches = 'tight')
        plt.show()
        

    """comparing betahat from RK4(5) to betahat from perturbed and analytical expression"""
    if fw_err == False and b_per is not None and b_analytical is not None:
        b_pert = b_per
        save_path = f"Plots/{base_name}_beta_vs_analytical_perturbed.png"
        plt.figure()
        plt.plot(t[r >= 0.1] , b[r >= 0.1] , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t[r >= 0.1] , b_pert[r >= 0.1] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\beta}$")
        plt.title(fr"$\hat{{\beta}}$ from RK4(5) and perturbed solution, {material}")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.show()

        plt.figure()
        plt.plot(t[r >= 0.1] , b[r >= 0.1] , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t[r >= 0.1] , b_analytical[r >= 0.1] , color = "orange" , linestyle = "--" , label = r"Analytical $\hat{\beta}$")
        plt.title(fr"$\hat{{\beta}}$ from RK4(5) and analytical solution , {material}")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.show()
        
        plt.figure()
        plt.plot(t[r >= 0.1] , b[r >= 0.1] , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t[r >= 0.1] , b_pert[r >= 0.1] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\beta}$")
        plt.plot(t[r >= 0.1] , b_analytical[r >= 0.1] , color = "orange" , linestyle = "--" , label = r"Analytical $\hat{\beta}$")
        plt.title(fr"$\hat{{\beta}}$ from RK4(5), perturbed and analytical solution, {material}")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
        plt.show()
    
    """comparing relative forward errors in beta hat from RK45-perturbed expression and 
    RK45-analytical expression"""
    if fw_err == True and b_per is not None and b_analytical is not None:
        b_pert = b_per
        save_path1 = f"Plots/{base_name}_beta_rel_fw_err.png"
        fw_err_RK45_per = np.abs(b - b_pert) / np.abs(b)
        fw_err_RK45_analytical = np.abs(b - b_analytical) / np.abs(b)
        plt.plot(t[r >= 0.1] , fw_err_RK45_per[r >= 0.1] , color = "blue" , label = "Rel fw error RK4(5) vs perturbed")
        plt.plot(t[r >= 0.1] , fw_err_RK45_analytical[r >= 0.1] , color = "red" , linestyle = "--" , label = "Rel fw error RK4(5) vs analytical")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"Relative forward error")
        plt.title(f"Relative forward error, RK4(5) vs perturbed and RK4(5) vs analytical, {material}")

        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')
        plt.show()
    
    if b_per is None and b_analytical is None and fw_err == False:
        # plt.plot(t[r >= 0.1] , b[r >= 0.1] , color = "blue")
        plt.plot(t , b , color = "blue")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.title(fr"$\hat{{\beta}}$ for {material}")
        plt.savefig(save_path2 , dpi = 300 , bbox_inches = 'tight')
        plt.show()

"""plots and compares energies between RK4(5) and Leapfrog solver"""
def energy_plot(solver1 , solver2 , particle_obj , fw_err = False):
    """input: solver1 (.npz), consisting of x1, y1 , vx1 , vx2 , m1 , b_vals1, RK4(5) solver
              t_arr (tuple), consisting of dt and t_tot, time at which solver1 and 2 have been 
                             evaluated
              solver2 (.npz), x2 , y2 , vx1 , vx2 , m2 , b_vals2 , Leapfrog solver

       returns: none"""
    
    res1 = np.load(solver1)
    x1 , y1 , vx1 , vy1 , m1 , b_vals1 , t = [res1[k] for k in ("x","y","vx","vy","m","b", "t")] #unpacking solver1

    res2 = np.load(solver2)
    x2 , y2 , vx2 , vy2 , m2 , b_vals2 , tl = [res2[k] for k in ("x","y","vx","vy","m","b", "t")] #unpacking solver2

    base_name = os.path.splitext(os.path.basename(solver1))[0]
    save_path = f"Plots/{base_name}_energy_solver_comps.png"

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

    if fw_err == False:
        plt.plot(0 , 0 , c = "black" , linestyle = "-" , label = "RK4(5)")
        plt.plot(0 , 0 , c = "black" , linestyle = "--" , label = "Leapfrog")
        plt.plot(t , kinetic1 , label = "Kinetic" , color = "blue" , linewidth = 2)
        plt.plot(t , potential1 , label = "Potential" , color = "orange" , linewidth = 2)
        plt.plot(t , tot1 , label = "Total" , color = "teal" , linewidth = 2)
        plt.plot(tl , kinetic2 , label = "Kinetic" , color = "red" , linestyle = "--" , linewidth = 2)
        plt.plot(tl , potential2 , label = "Potential" , color = "purple" , linestyle = "--" , linewidth = 2)
        plt.plot(tl , tot2 , label = "Total" , color = "pink" , linestyle = "--" , linewidth = 2)
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel("Energy")
        plt.title("RK4(5) vs leapfrog energies")
        plt.legend(loc = "upper right" ,
               bbox_to_anchor = (1.0 , 0.95))
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
        plt.show()
        
    if fw_err == True: #plots RK4(5) sols and relative forward errors between RK4(5) and Leapfrog
        save_path1 = f"Plots/{base_name}_energy_rel_fw_err.png"
        kinetic2 = np.interp(t , tl , kinetic2)
        potential2 = np.interp(t , tl , potential2)
        tot2 = np.interp(t , tl , tot2)

        err_kin = np.abs(kinetic1 - kinetic2) / np.abs(kinetic1)
        err_pot = np.abs(potential1 - potential2) / np.abs(potential1)
        err_tot = np.abs(tot1 - tot2) / np.abs(tot1)

        plt.figure()
        plt.plot(t , err_kin , label = "Kinetic" , color = "blue")
        plt.plot(t , err_pot  , label = "Potential" , color = "orange")
        # plt.plot(t , err_tot , label = "Total energy" , color = "teal" )
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel("Relative forward error")
        plt.title("Relative forward error between RK4(5) and leapfrog")
        plt.legend(loc = "upper right" ,
               bbox_to_anchor = (1.0 , 0.95))
        # plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')
        plt.show()
        
"""plots general beta curves for silicate and carbon"""
def beta_curves(interp = False , material = "silicate" , comp = False , scaled = False):
    """input: None
    
       returns: None"""
    
    """compare interpolated function with true curve for silicate values"""
    if interp == True: 
        size = sil_size if material == "silicate" else car_size
        beta = sil_betaval if material == "silicate" else car_betaval

        interp = pchip(size , beta)

        plt.plot(size , beta , linestyle = "-" , label = "True curve")
        plt.plot(size , interp(size) , linestyle = "--" , label = "Interpolated function")
        plt.xlabel(r"Particle size (m)")
        plt.ylabel(r"$\beta$")
        plt.title(fr"{material.capitalize()} $\beta$ versus particle size")
        plt.legend()
        plt.savefig(f"Plots/{material}_beta_interpolation_curve.png" , dpi = 300 , bbox_inches = 'tight')
        
    
    """comparing real silicate and carbon beta curve"""
    if comp == True:
        plt.xscale("log")
        plt.yscale("log")
        plt.plot(sil_size , sil_betaval , color = "red" , linestyle = "-" , label = "Silicate")
        plt.plot(car_size , car_betaval , color = "blue" , linestyle = "--" , label = "Carbon")
        plt.xlabel(r"Particle size (m)")
        plt.ylabel(r"$\beta$")
        plt.title(r"Silicate and carbon, $\beta$ versus particle size")
        plt.legend()
        plt.savefig(f"Plots/beta_interpolation_curve_silicate_carbon.png" , dpi = 300 , bbox_inches = 'tight')
    
    """example scaled curve silicate"""
    if scaled:
        ref_size = 10.33226 * 10**(-6)
        ref_B = 0.0163

        interp = pchip(sil_size , sil_betaval)
        plt.plot(sil_size / ref_size , interp(sil_size) / ref_B)
        plt.xlabel(r"Particle size (m)")
        plt.ylabel(r"$\hat{\beta}$")
        plt.title(r"Silicate $\hat{\beta}$ versus particle size")
        plt.savefig(f"Plots/silicate_betacurve_refscaling.png" , dpi = 300 , bbox_inches = 'tight')


    plt.show()

"""plots theoretical PR and sputtering lifetimes, silicate and carbon, all solar wind conditions"""
def PR_spu_lifetime():
    """input: None
    
    returns: None"""

    material = ["silicate" , "carbon"]
    sw_conds = ["slow" , "fast" , "CME"]
    tsp_vals = {m: {} for m in material}

    ls = {"silicate" : "-" , 
          "carbon" : "--"}
    
    cl = {"slow" : "green" ,
          "fast" : "red" ,
          "CME" : "blue"}
    
    for sw in sw_conds:
        par_sil = dust_properties("silicate" , sw , size = None , size_range = (sil_size , sil_betaval))
        t_sp_sil = par_sil.sputtering_lifetime()
        tsp_vals["silicate"][sw] = t_sp_sil

        par_carb = dust_properties("carbon" , sw , size = None , size_range = (car_size , car_betaval))
        t_sp_car = par_carb.sputtering_lifetime()
        tsp_vals["carbon"][sw] = t_sp_car

    for mat in material: 
        fig , ax = plt.subplots()
        size = sil_size if mat == "silicate" else car_size
        pr = tau_sil if mat == "silicate" else tau_car

        for sw, vals in tsp_vals[mat].items():
            ax.plot(size , vals,
                     color = cl[sw] , linestyle = ls[mat])
            
        ax.plot(0 , 0 , c = cl[sw])
        ax.plot(size , pr , color = "purple" , linestyle = ls[mat])

        purple_patch = mpatches.Patch(color = "purple" , label = "PR")
        blue_patch = mpatches.Patch(color = "blue" , label = "CME")
        red_patch = mpatches.Patch(color = "red" , label = "Fast")
        green_patch = mpatches.Patch(color = "green" , label = "Slow")

        handles, labels = plt.gca().get_legend_handles_labels()
        handles.extend([purple_patch , blue_patch , red_patch , green_patch])

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(0.1 , 10**8)

        ax.set_xlabel(r"Particle size (m)")
        ax.set_ylabel(r"Lifetime (years)")

        ax.set_title(f"{mat.capitalize()} PR and sputtering lifetimes")
        ax.legend(handles = handles)
        fig.savefig(f"Plots/{mat}_PR_sputtering_lifetime_theoretical.png", dpi = 300 , bbox_inches = 'tight')

    plt.show()
        
def PR_spu_lifetime_separate(lifetime_effects = "both"):
    material = ["silicate" , "carbon"]
    sw_conds = ["slow" , "fast" , "CME"]
    tsp_vals = {m: {} for m in material}

    ls = {"silicate" : "-" , 
          "carbon" : "--"}
    
    cl = {"slow" : "green" ,
          "fast" : "red" ,
          "CME" : "blue"}
    
    markers = {"pr" : "o" ,
               "sputtering" : "^" ,
               "both" : "x"}
    
    for sw in sw_conds:
        par_sil = dust_properties("silicate" , sw , size = None , size_range = (sil_size , sil_betaval))
        t_sp_sil = par_sil.sputtering_lifetime()
        tsp_vals["silicate"][sw] = t_sp_sil

        par_carb = dust_properties("carbon" , sw , size = None , size_range = (car_size , car_betaval))
        t_sp_car = par_carb.sputtering_lifetime()
        tsp_vals["carbon"][sw] = t_sp_car

    for mat in material: 
        fig , ax = plt.subplots()
        size = sil_size if mat == "silicate" else car_size
        
        pr = tau_sil if mat == "silicate" else tau_car

        if mat == "carbon":
                ax.axvspan(0.01516 * 10**(-6) , 0.54840 * 10**(-6) , color = "blue" , alpha = 0.1 , label = "B>1")

        for sw, vals in tsp_vals[mat].items():
            
            for key , value in true_lifetime_variableeps.items():

                ax.scatter(value["size"] , value[mat][lifetime_effects][sw] , c = cl[sw] , marker = markers[lifetime_effects])

            ax.plot(0 , 0 , c = cl[sw])

            ax.plot(size , vals,
                     color = cl[sw] , linestyle = ls[mat])
 
        ax.plot(size , pr , color = "purple" , linestyle = ls[mat])

        # ax.plot(0 , 0 , c = "black" , linestyle = ls[mat] , label = f"{mat.capitalize()}")
        ax.scatter(0 , 0 , c = "black" , marker = markers[lifetime_effects] , label = f"{lifetime_effects.capitalize()} numerical")

        purple_patch = mpatches.Patch(color = "purple" , label = "PR")
        blue_patch = mpatches.Patch(color = "blue" , label = "CME")
        red_patch = mpatches.Patch(color = "red" , label = "Fast")
        green_patch = mpatches.Patch(color = "green" , label = "Slow")

        handles, labels = plt.gca().get_legend_handles_labels()
        handles.extend([purple_patch , blue_patch , red_patch , green_patch])

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(0.1 , 10**8)

        ax.set_xlabel(r"Particle size (m)")
        ax.set_ylabel(r"Lifetime (years)")

       
        ax.set_title(f"{mat.capitalize()} PR and sputtering lifetimes theoretical vs {lifetime_effects} numerical" , pad = 20)
        ax.legend(handles = handles , fontsize = 8)
        fig.savefig(f"Plots/{mat}_PR_sputtering_lifetime_separate_{lifetime_effects}.png" , dpi = 300 , bbox_inches = 'tight')

    plt.show()
    

def v_theta(file_path , pert = None , material = None):

    res = np.load(file_path)
    x , y , vx , vy , _ , _ , t = [res[k] for k in ("x","y","vx","vy","m","b" , "t")] #unpacking file_path
    r = np.sqrt(x**2 + y**2)

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_vtheta.png"

    theta = np.atan2(y , x)
    theta = np.unwrap(theta)

    vtheta = -vx * np.sin(theta) + vy * np.cos(theta)

    if pert is not None:
        vtheta_pert = pert
        save_path1 = f"Plots/{base_name}_vtheta_pert_comps.png"

        plt.plot(t[r >= 0.1] , vtheta[r >= 0.1] , color = "blue" , label = "RK4(5)")
        plt.plot(t[r >= 0.1] , vtheta_pert[r >= 0.1] , color = "red" , linestyle = "--" , label = "Perturbed")


        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_{\theta}$")
        plt.title(fr"Perturbed and numerical $\hat{{v}}_{{\theta}}$, {material}")
        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')
    
    if pert is None:
        # plt.plot(t[r >= 0.1] , vtheta[r >= 0.1] , color = "blue")
        plt.plot(t , vtheta , color = "blue")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_{\theta}$")
        
        plt.title(fr"$\hat{{v}}_{{\theta}}$ for {material}")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    plt.show()

"""Plotting eccentricity from math and scaled expression"""
def ecc_math(file_path , pert = None):
    """input: x (array), x position
              y (array), y position
              vx (array), x velocity
              vy (array), y velocity
              beta (array), beta values
              t (array), time values
              type (string), default: None, "numerical" if using numerical solved params
              else "perturbed", if using perturbed params"""
    
    res = np.load(file_path)
    x , y , _ , _ , _ , _ , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b" , "t")] #unpacking file_path
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_ecc_ellipse_params.png"

    theta = np.atan2(y , x) #angle
    theta_cont = np.unwrap(theta) #continuous angle
    # orbit = np.floor(theta_cont / (2 * np.pi)).astype(int) #iterating over orbits

    plt.xlabel(r"$\hat{\theta}$ / $2\pi$")
    plt.ylabel(r"e value")
    plt.title("Eccentricity as function of orbits, elliptical parameters")

    if pert is not None:
        save_path1 = f"Plots/{base_name}_ecc_ellipse_params_pert_comps.png"
        ecc_num , ecc_pert , orb = ecc_calcs(x , y , pert) #mathematical eccentricity
        
        plt.plot(orb , ecc_pert , color = "red" , linestyle = "--" , label = "Perturbed")
        plt.plot(orb , ecc_num , color = "blue" , label = "Numerical")
        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')

    else:
        ecc_num , _ , orb = ecc_calcs(x , y) #mathematical eccentricity
        plt.plot(orb , ecc_num , color = "blue")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    plt.show()

def ecc_sc(file_path , B , pert = None):
    """input: file_path (str), containing file name, path and .filetype
              B (float), initial beta value corresponding to particle considered
              pert (tuple), rpert, omegapert, vrpert"""
    
    res = np.load(file_path)
    x , y , vx , vy , m , b , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b" , "t")] #unpacking file_path
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_ecc_scaled.png"

    theta = np.atan2(y , x) #angle
    theta_cont = np.unwrap(theta) #continuous angle
    # orbit = np.floor(theta_cont / (2 * np.pi)).astype(int) #iterating over orbits

    plt.title("Scaled eccentricity as function of orbits")
    plt.xlabel(r"$\hat{\theta}$ / $2\pi$")
    plt.ylabel(r"e value")

    if pert is not None:
        save_path1 = f"Plots/{base_name}_ecc_scaled_pert_comps.png"
        ecc_num , ecc_pert , orb = ecc_scaled((x , y , vx , vy , m , b , t) , B , pert) #mathematical eccentricity
        
        plt.plot(orb[:-1] , ecc_pert[:-1] , color = "red" , linestyle = "--" , label = "Perturbed")
        plt.plot(orb[:-1] , ecc_num[:-1] , color = "blue" , label = "Numerical")
        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')

    else:
        ecc_num , _ , orb = ecc_scaled((x , y , vx , vy , m , b , t) , B) #mathematical eccentricity
        plt.plot(orb[:-1] , ecc_num[:-1] , color = "blue")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    
    
    plt.show()

def mass_plot(file_path , file_path_comp , material):
    res = np.load(file_path)
    _ , _ , _ , _ , m , _ , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")] 
    base_name = os.path.splitext(os.path.basename(file_path))[0]

    res1 = np.load(file_path_comp)
    _ , _ , _ , _ , m1 , _ , t1 = [res1[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    base_name1 = os.path.splitext(os.path.basename(file_path_comp))[0]
    save_path1 = f"Plots/{base_name}_m_vs_{base_name1}_{material}.png"

    plt.plot(t , m , color = "blue", label = r"$\hat{m}$, constant $f_{sw}$")
    plt.plot(t1 , m1 , color = "red" , linestyle = "--" , label = r"$\hat{m}$, $f_{sw}(\hat{r})$")
    plt.title(fr"$\hat{{m}}$, constant $f_{{sw}}$ and $f_{{sw}}(\hat{{r}})$, {material}")
    plt.xlabel(r"$\hat{t}$")
    plt.ylabel(r"$\hat{m}$")
    plt.legend()

    plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')
        
    
    plt.show()

def eval_sizes():
    sizes = [1.54079 * 10**(-6) , 0.17508 * 10**(-6) , 0.04259 * 10**(-6)]
    sil_betas = [0.1235 , 0.8560 , 0.2098]
    car_betas = [0.2646 , 3.0589 , 1.6179]

    plt.plot(sil_size , sil_betaval , linestyle = "-" , color = "red" , label = "Silicate")
    plt.plot(car_size , car_betaval , linestyle = "--" , color = "blue" , label = "Carbon")

    plt.scatter(sizes , sil_betas , marker = "o" , color = "C4" , zorder = 4)
    plt.scatter(sizes , car_betas , marker = "o" , color = "C2" , zorder = 4)

    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Particle size (m)")
    plt.ylabel(r"$\beta$")
    plt.title(r"Initial particle sizes and $\beta$ values")
    plt.legend()
    plt.savefig("Plots/beta_referencesizes.png" , dpi = 300 , bbox_inches = 'tight')

    plt.show()

if __name__ == "__main__":
    
    filepath = "Files/rk45_t7_small_silicate_slowsw.npz"
    par = dust_properties("silicate" , "slow" , "small")
    res = np.load(filepath)
    x , y , _ , _ , m , b , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    
    p = perturbed_functions(par , t , b)
    
    rpert = p.rad()

    rhat_comps(filepath , "silicate" , None , rpert)
    
    
    
    


    

    
    
    
    
    