import matplotlib.pyplot as plt
import numpy as np
from dust_properties import dust_properties
from energy import tot_energy
from config import (car_betaval_bound , car_size_bound , init_vals , sil_beta , car_beta , t6 , t7 ,
                    t8 , t9 , t10 , sil_size , sil_betaval , car_size , car_betaval , sil_PR , car_PR 
                    , tau_car , tau_sil , size_to_mass)
from forces_scaled import betahat
from scipy.interpolate import PchipInterpolator as pchip
import os
# from lifetime_calcs import true_lifetime , true_lifetime_variableeps
from new_num_lifetimes import true_lifetime , true_lifetime_variableeps
import matplotlib.patches as mpatches
from eccentricity import ecc_calcs , ecc_scaled
from new_pert_lifetimes import pert_lifetime
import matplotlib
from pert_variable_eps import PerturbationSolver

"""plotting params to adjust font sizes"""
plt.rcParams.update({"font.size" : 14 ,
    "axes.labelsize" : 14 ,
    "axes.titlesize" : 14 ,
    "xtick.labelsize" : 12 ,
    "ytick.labelsize" : 12 ,
    "legend.fontsize" : 12})

"""Epsilon for all possible B values, carbon and silicate"""
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
        fg , ax = plt.subplots(figsize = (5 , 4))
        for sw_cond , styles in sw.items():
            
            size_ranges = (sil_size , sil_betaval) if mat == "silicate" else (car_size_bound , car_betaval_bound)
            par = dust_properties(mat , sw_cond , init_dist = 1.0 , size = None , size_range = size_ranges)
            init_beta = par.B
            epsilon = par.eps()
            delta = par.delta
            
            ax.plot(init_beta , epsilon , color = styles["color"] , linestyle = "-")
            ax.plot(init_beta , delta , color = "purple" , linestyle = "-")
            
            #Arrow code from Dietrich in Stackoverflow
            ax.annotate('', xy=(1.03 , -0.03), xycoords='axes fraction', xytext=(1.03 , 1), 
                    arrowprops=dict(arrowstyle="->", color='black'))
            ax.text(1.05, 0.5,r"${a_0}$",transform=ax.transAxes,va='center')
            fg.canvas.draw()
            
        purple_patch = mpatches.Patch(color = "purple" , label = r"${\delta}$")
        blue_patch = mpatches.Patch(color = "blue" , label = r"${\epsilon}$ CME")
        red_patch = mpatches.Patch(color = "red" , label = r"${\epsilon}$ Fast")
        green_patch = mpatches.Patch(color = "green" , label = r"${\epsilon}$ Slow")

        handles, labels = plt.gca().get_legend_handles_labels()
        handles.extend([purple_patch , blue_patch , red_patch , green_patch])    
        plt.loglog()        
        plt.title(fr"{mat.capitalize()}, ${{\epsilon}}$ and ${{\delta}}$ vs B")
        plt.xlabel("B")
        plt.ylabel("Value")
        plt.ylim(10**(-9) , 1)
        plt.legend(handles = handles)
        fg.subplots_adjust(right=0.9)
        plt.savefig(f"Plots/{mat}_epsilonvsbeta.png", dpi = 300 , bbox_inches = 'tight')
        plt.show()
    
"""thetahat RK4(5) only or RK4(5) vs perturbed"""
def thetahat_comps(file_path , pert = None , material = None):
    """input: file_path (string), path for RK4(5) simulation file
              pert (array), default:None, else array containing perturbed values
              material (string), default:None, silicate or carbon
              
        returns: none   """
    
    res = np.load(file_path)
    x , y , _ , _ , _ , _ , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b" , "t")]
    r = np.sqrt(x**2 + y**2)

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_theta.png"

    theta1 = np.atan2(y , x) #rk45 thetahat
    theta1 = np.unwrap(theta1) #avoiding discontinuities in thetahat

    """comparing thetahat RK4(5) and perturbed"""
    if pert is not None:
        theta_per = pert
        save_path2 = f"Plots/{base_name}_theta_vs_perturbed.png"
        theta_per = np.unwrap(theta_per) #removing discontinuities

        plt.plot(t[r >= 0.1] , theta1[r >= 0.1] , color = "blue" , label = r"RK4(5) $\hat{\theta}$")
        plt.plot(t[r >= 0.1] , theta_per[r >= 0.1] , color = "red" , linestyle = "--" 
                 , label = r"Perturbed $\hat{\theta}$")
        
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\theta}$")
        plt.title(rf"$\hat{{\theta}}$ from RK4(5) and perturbed solution, {material}")
        plt.legend()
        plt.savefig(save_path2 , dpi = 300 , bbox_inches = 'tight')
    
    """only RK4(5)"""
    if pert is None:
        plt.plot(t[r >= 0.1] , theta1[r >= 0.1] , color = "blue")
        # plt.plot(t , theta1 , color = "blue")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\theta}$")

        plt.title(rf"$\hat{{\theta}}$ for {material}")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    plt.show()

"""rhat RK4(5) only or comparison RK4(5) with perturbed"""
def rhat_comps(file_path , pert = None , material = None , mat_comp = None):
    """input: file_path (string), path for RK4(5) simulation file
              pert (array), default:None, else array containing perturbed values
              material (string), default:None, silicate or carbon
              
        returns: none   """
    
    res = np.load(file_path)
    x , y , _ , _ , _ , _ , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    r = np.sqrt(x**2 + y**2) #r hat
    
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_r.png"

    """RK4(5) vs perturbed"""
    if pert is not None: #comparing RK4(5) with perturbed rhat
        plt.figure(figsize = (5 , 4))
        r_per,time = pert
        
        mask = round(len(r) / len(r_per))
        # r_per = r_per
        # t = t[::mask]
        # r = r[::mask]
        
        # rel_fw_err = np.abs(r - r_per) / np.abs(r)
        save_path2 = f"Plots/{base_name}_r_vs_perturbed.png"
        
        plt.plot(t , r , color = "blue" , label = r"RK4(5) $\hat{r}$")
        plt.plot(time , r_per , color = "red" , linestyle = "--" 
                 , label = r"Perturbed $\hat{r}$")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{r}$")
        plt.title(fr"{material.capitalize()} $\hat{{r}}$, RK4(5) and perturbed solution")
        plt.legend()
        plt.savefig(save_path2 , dpi = 300 , bbox_inches = 'tight')
    
    """only RK4(5)"""
    if pert is None and mat_comp is None:
        plt.figure(figsize = (5 , 4))
        plt.plot(t[r >= 0.1] , r[r >= 0.1] , color = "blue" , label = fr"$\hat{{r}}$ for {material}")
        # plt.plot(t , r , color = "blue" , label = fr"$\hat{{r}}$ for {material}")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{r}$")
        plt.title(fr"$\hat{{r}}$ for {material}")
        
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
    
    if mat_comp is not None:
        save_path1 = f"Plots/{base_name}_r_both_materials.png"
        file1 = mat_comp
        res1 = np.load(file1)
        x1 , y1 , _ , _ , _ , _ , t1 = [res1[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
        r1 = np.sqrt(x1**2 + y1**2)
        
        # plt.figure(figsize = (5 , 4))
        plt.plot(t , r , color = "blue" , label = r"Silicate $\hat{r}$")
        plt.plot(t1 , r1 , color = "red" , linestyle = "--" , label = r"Carbon $\hat{r}$")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{r}$")
        plt.title(r"$\hat{r}$ for silicate and carbon")
        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')
    
    plt.show()

""""vrhat RK4(5) only or RK4(5) vs perturbed"""
def vhat_comps(file_path , pert = None , material = None):
    """input: file_path (string), path for RK4(5) simulation file
              pert (array), default:None, else array containing perturbed values
              material (string), default:None, silicate or carbon
              
        returns: none   """

    res = np.load(file_path)
    x , y , vx , vy , _ , _ , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b" , "t")]
    r = np.sqrt(x**2 + y**2)

    theta_num = np.atan2(y , x) #thetahat
    theta_num = np.unwrap(theta_num) #avoiding discontinuities

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_vr.png"

    v_r = vx * np.cos(theta_num) + vy * np.sin(theta_num) #cartesian to radial vel
    
    """RK4(5) vs perturbed"""
    if pert is not None:
        v_per = pert
        
        save_path1 = f"Plots/{base_name}_vr_vs_perturbed.png"
        plt.plot(t[r >= 0.1] , v_r[r >= 0.1] , color = "blue" , label = r"RK4(5) $\hat{v}_r$")
        plt.plot(t[r >= 0.1] , v_per[r >= 0.1] , color = "red" , linestyle = "--" 
                 , label = r"Perturbed $\hat{v}_r$")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_r$")

        plt.legend()
        plt.title(rf"{material.capitalize()} $\hat{{v}}_r$, RK4(5) and perturbed solution")
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')
    
    """only RK4(5)"""
    if pert is None:
        # plt.plot(t[r >= 0.1] , v_r[r >= 0.1] , color = "blue")
        plt.plot(t , v_r , color = "blue")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_r$")
        plt.title(fr"$\hat{{v}}_r$ for {material}")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    plt.show()

"""omegahat RK4(5) only or RK4(5) vs perturbed"""
def omegahat_comps(file_path , pert = None , material = None , mat_comp = None):
    """input: file_path (string), path for RK4(5) simulation file
              pert (array), default:None, else array containing perturbed values
              material (string), default:None, silicate or carbon
              
        returns: none   """
    
    res = np.load(file_path)
    x , y , vx , vy , _ , _ , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    r = np.sqrt(x**2 + y**2)

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_omega.png"
    
    theta_num = np.atan2(y , x) #thetahat RK4(5)
    theta_num = np.unwrap(theta_num) #avoiding discontinuities

    angvel_num = (-vx * np.sin(theta_num) + vy * np.cos(theta_num)) / r    

    """RK4(5) vs perturbed"""
    if pert is not None:
        plt.figure(figsize = (5 , 4))
        angvel, time = pert

        # t = t[:len(angvel)]
        # angvel_num = angvel_num[:len(angvel)]
        # r = r[:len(angvel)]

        plt.plot(t , angvel_num , color = "blue" , label = r"RK4(5) $\hat{\omega}$")
        plt.plot(time , angvel , color = "red" , linestyle = "--" 
                 , label = r"Perturbed $\hat{\omega}$")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\omega}$")
        plt.title(rf"{material.capitalize()} $\hat{{\omega}}$, RK4(5) and perturbed solution")

        plt.legend()
        save_path = f"Plots/{base_name}_omega_vs_perturbed.png"
        
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    """only RK4(5)"""
    if pert is None and mat_comp is None:
        plt.figure(figsize = (5 , 4))
        plt.plot(t[r >= 0.1] , angvel_num[r >= 0.1] , color = "blue")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\omega}$")
        
        plt.title(rf"$\hat{{\omega}}$ for {material}")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
    
    if mat_comp is not None:
        save_path1 = f"Plots/{base_name}_omega_both_materials.png"
        file1 = mat_comp
        res1 = np.load(file1)
        x1 , y1 , vx1 , vy1 , _ , _ , t1 = [res1[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
        r1 = np.sqrt(x1**2 + y1**2)
        theta_num1 = np.atan2(y1 , x1) #thetahat RK4(5)
        theta_num1 = np.unwrap(theta_num1) #avoiding discontinuities

        angvel_num1 = (-vx1 * np.sin(theta_num1) + vy1 * np.cos(theta_num1)) / r1   
        
        # plt.figure(figsize = (5 , 4))
        plt.plot(t , angvel_num , color = "blue" , label = r"Silicate $\hat{\omega}$")
        plt.plot(t1 , angvel_num1 , color = "red" , linestyle = "--" , label = r"Carbon $\hat{\omega}$")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\omega}$")
        plt.title(r"$\hat{\omega}$ for silicate and carbon")
        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')

    plt.show()

"""betahat values RK4(5) only or RK4(5) vs perturbed"""
def b_plot(file_path , b_per = None , material = None , mat_comp = None):
    """input: file_path (string), path for RK4(5) simulation file
              b_per (array), default:None, else array containing perturbed values
              material (string), default:None, silicate or carbon
              
        returns: none   """
    
    res = np.load(file_path)
    x , y , _ , _ , _ , b , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    r = np.sqrt(x**2 + y**2)

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path2 = f"Plots/{base_name}_beta.png"

    """RK4(5) vs perturbed"""
    if b_per is not None:
        plt.figure(figsize = (5 , 4))

        b_pert , time = b_per
        # t = t[:len(b_pert)]
        # b = b[:len(b_pert)]
        # r = r[:len(b_pert)]

        save_path = f"Plots/{base_name}_beta_vs_perturbed.png"

        plt.plot(t , b , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(time , b_pert , color = "red" , linestyle = "--" 
                 , label = r"Perturbed $\hat{\beta}$")
        plt.title(fr"{material.capitalize()} $\hat{{\beta}}$, RK4(5) and perturbed solution")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
        plt.show()
    
    """only RK4(5)"""
    if b_per is None and mat_comp is None:
        plt.figure(figsize = (5 , 4))
        plt.plot(t[r >= 0.1] , b[r >= 0.1] , color = "blue")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.title(fr"$\hat{{\beta}}$ for {material}")

        plt.savefig(save_path2 , dpi = 300 , bbox_inches = 'tight')
        plt.show()
    
    if mat_comp is not None:
        save_path1 = f"Plots/{base_name}_beta_both_materials.png"
        file1 = mat_comp
        res1 = np.load(file1)
        x1 , y1 , vx1 , vy1 , _ , b1 , t1 = [res1[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
        r1 = np.sqrt(x1**2 + y1**2) 
        
        # plt.figure(figsize = (5 , 4))
        plt.plot(t , b , color = "blue" , label = r"Silicate $\hat{\beta}$")
        plt.plot(t1 , b1 , color = "red" , linestyle = "--" , label = r"Carbon $\hat{\beta}$")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ for silicate and carbon")
        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')
        plt.show()

"""plots and compares energies between RK4(5) and Leapfrog solver"""
def energy_plot(solver1 , solver2 , particle_obj , fw_err = False):
    """input: solver1 (string), consisting of x1, y1 , vx1 , vx2 , m1 , b_vals1, t1, RK4(5) solver
              solver2 (string), x2 , y2 , vx1 , vx2 , m2 , b_vals2 , t2 Leapfrog solver
              particle_obj (instance), containing current particle information
              fw_err (bool), default:False, choose if calculating relative forward error

       returns: none"""
    
    res1 = np.load(solver1)
    x1 , y1 , vx1 , vy1 , m1 , b_vals1 , t1 = [res1[k] for k in ("x","y","vx","vy","m","b", "t")] #unpacking solver1

    res2 = np.load(solver2)
    x2 , y2 , vx2 , vy2 , m2 , b_vals2 , t2 = [res2[k] for k in ("x","y","vx","vy","m","b", "t")] #unpacking solver2

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
        plt.plot(t1 , kinetic1 , label = "Kinetic" , color = "blue" , linewidth = 2)
        plt.plot(t1 , potential1 , label = "Potential" , color = "orange" , linewidth = 2)
        plt.plot(t1 , tot1 , label = "Total" , color = "teal" , linewidth = 2)
        plt.plot(t2 , kinetic2 , label = "Kinetic" , color = "red" , linestyle = "--" , linewidth = 2)
        plt.plot(t2 , potential2 , label = "Potential" , color = "purple" , linestyle = "--" 
                 , linewidth = 2)
        plt.plot(t2 , tot2 , label = "Total" , color = "pink" , linestyle = "--" , linewidth = 2)
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel("Energy")
        plt.title("RK4(5) vs leapfrog energies")
        plt.legend(loc = "upper right" ,
               bbox_to_anchor = (1.0 , 0.95))
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
        plt.show()
        
    if fw_err == True: #plots RK4(5) sols and relative forward errors between RK4(5) and Leapfrog
        kinetic2 = np.interp(t , t2 , kinetic2)
        potential2 = np.interp(t , t2 , potential2)
        tot2 = np.interp(t , t2 , tot2)

        err_kin = np.abs(kinetic1 - kinetic2) / np.abs(kinetic1)
        err_pot = np.abs(potential1 - potential2) / np.abs(potential1)
        err_tot = np.abs(tot1 - tot2) / np.abs(tot1)
        
"""plots general beta curves for silicate and carbon"""
def beta_curves(interp = False , material = "silicate" , comp = False , pert = False):
    """input: interp (bool), default:False, choose if wanting to compare an interpolated function
                            with the experimental curve
              material (string), silicate or carbon material to consider
              comp (bool), default:False, choose if wanting only experimental curves
    
       returns: None"""
    
    """compare interpolated function with experimental curve for material values"""
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
        rel_fw_err = np.abs(beta - interp(size)) / np.abs(beta)
        plt.savefig(f"Plots/{material}_beta_interpolation_curve.png" , dpi = 300 , bbox_inches = 'tight')
        
    
    """experimental silicate and carbon beta curve"""
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

    if pert == True:
        ref_size = 50 * 10**(-6)
        bpert_sil = (sil_size / ref_size)**(-1)
        bpert_car = (car_size / ref_size)**(-1)
        refB_sil = 0.0032
        refB_car = 0.0063
        plt.xscale("log")
        plt.yscale("log")
        plt.title(r"Silicate and carbon, numerical and analytical $\beta$")
        plt.plot(sil_size , sil_betaval , color = "red" , linestyle = "-")
        plt.plot(sil_size , bpert_sil * refB_sil , color = "red" , linestyle = "--")
        
        plt.plot(car_size , car_betaval , color = "blue" , linestyle = "-")
        plt.plot(car_size  , bpert_car * refB_car , color = "blue" , linestyle = "--")
        plt.xlabel(r"Particle size (m)")
        plt.ylabel(r"$\beta$")
        plt.plot(0 , 0 , c = "black" , linestyle = "--" , label = r"Analytical $\beta$")
        plt.plot(0 , 0 , c = "black" , linestyle = "-" , label = r"Experimental $\beta$")
        red_patch = mpatches.Patch(color = "red" , label = "Silicate")
        blue_patch = mpatches.Patch(color = "blue" , label = "Carbon")


        handles, labels = plt.gca().get_legend_handles_labels()
        handles.extend([red_patch , blue_patch])

        plt.legend(handles = handles , fontsize = 8)
        plt.savefig(f"Plots/beta_car_exp_analytical.png" , dpi = 300 , bbox_inches = 'tight')
        
        
    plt.show()

"""plots theoretical PR and sputtering lifetimes, silicate and carbon, all solar wind conditions"""
def PR_spu_lifetime_theo():
    """input: None
    
    returns: None"""

    material = ["silicate" , "carbon"]
    sw_conds = ["slow" , "fast" , "CME"]
    tsp_vals = {m: {} for m in material}
    
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

            ax.plot(0 , 0 , c = cl[sw])

            ax.plot(size , vals,
                     color = cl[sw] , linestyle = "-")
             
        ax.plot(size , pr , color = "purple" , linestyle = "-")

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

       
        ax.set_title(f"{mat.capitalize()} PR and sputtering lifetimes theoretical and numerical" 
                     , pad = 20)
        ax.legend(handles = handles , fontsize = 8)

        fig.savefig(f"Plots/{mat}__lifetimes_theoretical.png" , dpi = 300 
                        , bbox_inches = 'tight')

    plt.show()

"""plots theoretical PR and sputtering lifetimes, silicate and carbon, all solar wind conditions"""
def PR_spu_lifetime_num(file = true_lifetime , lifetime_effects = "both"):
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
                ax.axvspan(0.01516 * 10**(-6) , 0.54840 * 10**(-6) , color = "blue" , alpha = 0.1 
                           , label = "B>1")

        for sw, vals in tsp_vals[mat].items():
            
            for key , value in file.items():

                ax.scatter(value["size"] , value[mat][lifetime_effects][sw] , c = cl[sw] 
                           , marker = markers[lifetime_effects])

            ax.plot(0 , 0 , c = cl[sw])

            ax.plot(size , vals,
                     color = cl[sw] , linestyle = "-")
             
        ax.plot(size , pr , color = "purple" , linestyle = "-")
        
        ax.scatter(0 , 0 , c = "black" , marker = markers[lifetime_effects] 
                   , label = f"{lifetime_effects.capitalize()} effects")

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

       
        ax.set_title(f"{mat.capitalize()} PR and sputtering lifetimes theoretical and numerical" 
                     , pad = 20)
        ax.legend(handles = handles , fontsize = 8)

        fig.savefig(f"Plots/{mat}_{lifetime_effects}_lifetimes_cst_eps.png" , dpi = 300 
                        , bbox_inches = 'tight')

    plt.show()

"""Numerical PR and sputtering lifetimes vs pert"""    
def PR_spu_lifetime_pert(file = true_lifetime_variableeps , file_pert = pert_lifetime):
    """input: file (dictionary), default:true_lifetime_variableeps, other option true_lifetime
                                 choose if wanting constant or varying epsilon
              lifetime_effects (string), default: both, options pr, sputtering, which effects to consider
    
    returns: None"""

    material = ["silicate" , "carbon"]
    sw_conds = ["slow" , "fast" , "CME"]
    tsp_vals = {m: {} for m in material}
    
    cl = {"slow" : "green" ,
          "fast" : "red" ,
          "CME" : "blue"}
    
    markers = {"pr" : "o" ,
               "sputtering" : "^" ,
               "both" : "x"}
    pert_marker = matplotlib.markers.MarkerStyle('o' , fillstyle = 'none')
    
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
                ax.axvspan(0.01516 * 10**(-6) , 0.54840 * 10**(-6) , color = "blue" , alpha = 0.1 
                           , label = "B>1")

        for sw, vals in tsp_vals[mat].items():
            
            for key , value in file_pert.items():

                ax.scatter(value["size"] , value[mat]["both"][sw] , c = cl[sw] 
                           , marker = pert_marker)
                
            ax.plot(0 , 0 , c = cl[sw])

            ax.plot(size , vals,
                     color = cl[sw] , linestyle = "-")
            
            
            # for key , value in file.items():

            #     ax.scatter(value["size"] , value[mat][lifetime_effects][sw] , c = cl[sw] 
            #                , marker = markers[lifetime_effects])

            # ax.plot(0 , 0 , c = cl[sw])

            # ax.plot(size , vals,
            #          color = cl[sw] , linestyle = "-")
             
        ax.plot(size , pr , color = "purple" , linestyle = "-")
        
        # ax.scatter(0 , 0 , c = "black" , marker = markers[lifetime_effects] 
        #            , label = f"{lifetime_effects.capitalize()} effects")
        

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

       
        ax.set_title(f"{mat.capitalize()} PR and sputtering lifetimes theoretical and numerical" 
                     , pad = 20)
        
        ax.legend(handles = handles , fontsize = 8)

        fig.savefig(f"Plots/{mat}_botheffects_lifetimes_variable_eps_num_pert.png" , dpi = 300 
                        , bbox_inches = 'tight')


    plt.show()
    
"""vtheta RK4(5) or RK4(5) comparison with perturbed"""
def v_theta(file_path , pert = None , material = None):
    """input: file_path (string), path for RK4(5) simulation file
              pert (array), default:None, else array containing perturbed values
              material (string), default:None, silicate or carbon
              
        returns: none   """

    res = np.load(file_path)
    x , y , vx , vy , _ , _ , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b" , "t")] #unpacking file_path
    r = np.sqrt(x**2 + y**2)

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_vtheta.png"

    theta = np.atan2(y , x)
    theta = np.unwrap(theta)

    vtheta = -vx * np.sin(theta) + vy * np.cos(theta)

    """RK4(5) vs perturbed"""
    if pert is not None:
        vtheta_pert = pert
        save_path1 = f"Plots/{base_name}_vtheta_pert_comps.png"

        plt.plot(t[r >= 0.1] , vtheta[r >= 0.1] , color = "blue" , label = "RK4(5)")
        plt.plot(t[r >= 0.1] , vtheta_pert[r >= 0.1] , color = "red" , linestyle = "--" 
                 , label = "Perturbed")

        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_{\theta}$")
        plt.title(fr"Perturbed and numerical $\hat{{v}}_{{\theta}}$, {material}")
        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')
    
    """only RK4(5)"""
    if pert is None:
        plt.plot(t[r >= 0.1] , vtheta[r >= 0.1] , color = "blue")
        # plt.plot(t , vtheta , color = "blue")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_{\theta}$")
        
        plt.title(fr"$\hat{{v}}_{{\theta}}$ for {material}")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    plt.show()

"""Plotting eccentricity from elliptical definition, RK4(5) and pert"""
def ecc_math(file_path , pert = None):
    """input: file_path (string), path for RK4(5) simulation file
              pert (array), default:None, else array containing perturbed values

        returns: none   """
    
    res = np.load(file_path)
    x , y , _ , _ , _ , _ , _ = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b" , "t")] #unpacking file_path
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_ecc_ellipse_params.png"

    theta = np.atan2(y , x) #angle
    theta_cont = np.unwrap(theta) #continuous angle

    plt.xlabel(r"$\hat{\theta}$ / $2\pi$")
    plt.ylabel(r"e value")
    plt.title("Eccentricity as function of orbits, elliptical parameters")

    """RK4(5) vs perturbed"""
    if pert is not None:
        save_path1 = f"Plots/{base_name}_ecc_ellipse_params_pert_comps.png"
        ecc_num , ecc_pert , orb = ecc_calcs(x , y , pert) #ellipse def eccentricity
        
        plt.plot(orb , ecc_pert , color = "red" , linestyle = "--" , label = "Perturbed")
        plt.plot(orb , ecc_num , color = "blue" , label = "Numerical")
        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')

    """only RK4(5)"""
    if pert is None:
        ecc_num , _ , orb = ecc_calcs(x , y) #ellipse def eccentricity
        plt.plot(orb , ecc_num , color = "blue")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    plt.show()

"""Scaled eccentricity, RK4(5) and pert"""
def ecc_sc(file_path , particle_obj , pert = None):
    """input: file_path (str), containing file name, path and .filetype
              particle_obj (instance), information on dust properties
              pert (array), default:None, perturbed values 
        
        returns: None"""
    
    res = np.load(file_path)
    x , y , vx , vy , m , b , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b" , "t")] #unpacking file_path
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_ecc_scaled.png"

    theta = np.atan2(y , x) #angle
    theta_cont = np.unwrap(theta) #continuous angle

    plt.title("Scaled eccentricity as function of orbits")
    plt.xlabel(r"$\hat{\theta}$ / $2\pi$")
    plt.ylabel(r"e value")

    """RK4(5) vs perturbed"""
    if pert is not None:
        save_path1 = f"Plots/{base_name}_ecc_scaled_pert_comps.png"
        ecc_num , ecc_pert , orb = ecc_scaled((x , y , vx , vy , m , b , t) , particle_obj , pert) #mathematical eccentricity
        
        plt.plot(orb[:-1] , ecc_pert[:-1] , color = "red" , linestyle = "--" , label = "Perturbed")
        plt.plot(orb[:-1] , ecc_num[:-1] , color = "blue" , label = "Numerical")
        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')

    """only RK4(5)"""
    if pert is None:
        ecc_num , _ , orb = ecc_scaled((x , y , vx , vy , m , b , t) , particle_obj) #mathematical eccentricity
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
    
    plt.plot(sil_size , sil_betaval , linestyle = "-" , color = "blue" , label = "Silicate")
    plt.plot(car_size , car_betaval , linestyle = "--" , color = "red" , label = "Carbon")

    for label , i in init_vals.items():
        sizes = i["r"]
        sil_betas = i["B"]["silicate"]
        car_betas = i["B"]["carbon"]

        plt.scatter(sizes , sil_betas , marker = "o" , color = "C2" , zorder = 4)
        plt.scatter(sizes , car_betas , marker = "o" , color = "C4" , zorder = 4)

        # plt.annotate(label,(sizes, sil_betas),xytext=(3, 3),textcoords="offset points",fontsize=10)

        plt.annotate(label,(sizes, car_betas),xytext=(3, 3),textcoords="offset points",fontsize=10)

    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Particle size (m)")
    plt.ylabel(r"$\beta$")
    plt.title(r"Initial particle sizes and $\beta$ values")
    plt.legend()
    plt.savefig("Plots/beta_referencesizes.png" , dpi = 300 , bbox_inches = 'tight')

    plt.show()

if __name__ == "__main__":
    PR_spu_lifetime_pert()
    
    
    