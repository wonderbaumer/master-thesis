import matplotlib.pyplot as plt
import numpy as np
from forces import beta
from dust_properties import dust_properties
from pert_functions import perturbed_functions
from energy import tot_energy
from config import car_size_bound , car_betaval_bound , m_s , mA_S , sil_beta , car_beta , t5 , t6 , t7 , R , sil_size , sil_betaval , car_size , car_betaval , sil_PR , car_PR , sil_mass , car_mass
from forces_scaled import betahat
from scipy.interpolate import PchipInterpolator as pchip
from polar_to_cart import polar_to_cartesian
import os
from scipy.constants import G

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
    
    material = {"silicate" : "-" , 
                "carbon" : "--"}
    
    material = {"silicate" : "-" , 
                "carbon" : "--"}
    
    sw_conds = ["slow" , "fast" , "CME"]
    eps_vals = []
    delta_vals = []

    for key , value in material.items():
        par = dust_properties(key , "slow" , "all" , None)
        delta_vals.append({"material" : key ,
                           "delta": par.delta ,
                           "line" : value ,
                           "B" : par.B})
        
        for sw in sw_conds:
            par = dust_properties(key , sw , "all" , None)
            eps = par.epsilon
            drag_cst = par.K
            eps_vals.append({"material" : key ,
                             "sw_cond" : sw ,
                             "eps" : eps ,
                             "B" : par.B ,
                             "line" : value , 
                             "K" : drag_cst ,
                             "delta" : par.delta})   
        
    plt.figure()

    for item in eps_vals:
        m = item["material"]
        sw = item["sw_cond"]
        epsilon = item["eps"]
        Bval = item["B"]
        linestyle = item["line"]
        drag_cst = item["K"]
        dragtot = item["delta"]

        mask = np.abs((dragtot - epsilon / 5)) < 0.0000009

        print(m,sw)
        print(f"B {Bval[mask]},K {drag_cst[mask]}, epsilon {epsilon[mask]} , delta {dragtot[mask]}")
    #     plt.plot(Bval , epsilon , linestyle = linestyle , label=f"{m}, {sw} sw")
    # plt.xlabel(r"$B$")
    # plt.ylabel(r"${\epsilon}$")
    # plt.yscale("log")
    # plt.title(r"${\epsilon}$ vs ${B}$, corresponding to size range $1~\mathrm{nm} \text{–} 50~\mu\mathrm{m}$, silicate and carbon")
    # plt.legend(loc = "lower right" , prop = {"size" : 9} , borderpad = 0.05 , labelspacing = 0.2 , handlelength = 1.5)
    # plt.savefig("Plots/eps_vs_beta.png" , dpi = 300 , bbox_inches = 'tight')

    # plt.figure()
    # for item in delta_vals:
    #     m = item["material"]
    #     drag_term = item["delta"]
    #     Bval = item["B"]
    #     linestyle = item["line"]

    #     plt.plot(Bval , drag_term , linestyle = linestyle , label = f"{m}")
    # plt.xlabel(r"$B$")
    # plt.ylabel(r"${\delta}$")
    # plt.yscale("log")
    # plt.title(r"${\delta}$ vs ${B}$, corresponding to size range $1~\mathrm{nm} \text{–} 50~\mu\mathrm{m}$, silicate and carbon")
    # plt.legend(loc = "lower left")
    # plt.savefig("Plots/delta_vs_beta.png" , dpi = 300 , bbox_inches = 'tight')

    # plt.figure()
    # for i in eps_vals:
    #     Bval = i["B"]
    #     m = i["material"]
    #     epsilon = i["eps"]
    #     sw = i["sw_cond"]
    #     linestyle = i["line"]

    #     plt.plot(Bval , epsilon , linestyle = linestyle , label = f"{m}, {sw} sw eps")

    # for j in delta_vals:
    #     drag_term = j["delta"]
    #     Bval = j["B"]
    #     linestyle = j["line"]
    #     m = j["material"]

    #     plt.plot(Bval , drag_term , linestyle = linestyle , label = f"{m}, {sw} sw delta")
        
    # plt.xlabel(r"$B$")
    # plt.ylabel("Magnitude")
    # plt.yscale("log")
    # plt.title(r"${\delta}$ and ${\epsilon}$ comparison, silicate and carbon")
    # plt.legend(loc = "lower right" , prop = {"size" : 9} , borderpad = 0.05 , labelspacing = 0.2 , handlelength = 1.5)
    # plt.savefig("Plots/delta_vs_epsilon.png" , dpi = 300 , bbox_inches = 'tight')
    
    # plt.show()
    
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
    ymin , ymax = np.where(r <= 0.05) , np.where(r <= 1)

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
        
        plt.plot(t[::10] , theta1[::10] , color = "blue", label = r"RK4(5) $\hat{\theta}$")
        plt.plot(t[::10] , theta2[::10] , color = "red" , linestyle = "--" , label = r"Leapfrog $\hat{\theta}$")
        plt.title(r"$\hat{\theta}$ from RK4(5) and Leapfrog solution")
        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')

    """comparing thetahat from RK4(5) perturbed thetahat"""
    if pert is not None:
        theta_per = pert
        save_path2 = f"Plots/{base_name}_theta_vs_perturbed.png"
        theta_per = np.unwrap(theta_per) #removing discontinuities

        plt.plot(t[::10] , theta1[::10] , color = "blue" , label = r"RK4(5) $\hat{\theta}$")
        plt.plot(t[::10] , theta_per[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\theta}$")

        # if theta1[ymin][0] < theta_per[ymin][0]:
        #     plt.ylim(theta1[ymin][0] , theta1[ymax][-1])
        
        # else:
        #     plt.ylim(theta_per[ymin][0] , theta_per[ymax][-1])
        
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\theta}$")
        plt.title(r"$\hat{\theta}$ from RK4(5) and perturbed solution")
        plt.legend()
        # plt.savefig(save_path2 , dpi = 300 , bbox_inches = 'tight')
    
    if material == "silicate":
        plt.plot(t[::10] , theta1[::10])
        plt.ylim(theta1[ymin][0] , theta1[ymax][-1])
        plt.title(r"$\hat{\theta}$ for silicate, real $\hat{\beta}$")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
    
    if material == "carbon":
        plt.plot(t[::10] , theta1[::10])
        plt.title(r"$\hat{\theta}$ for carbon, real $\hat{\beta}$")
        plt.ylim(theta1[ymin][0] , theta1[ymax][-1])
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    plt.xlabel(r"$\hat{t}$")
    plt.ylabel(r"$\hat{\theta}$")
    plt.show()

"""plotting rhat as function of orbits, comparison of Leapfrog and RK4(5) solver, or RK4(5)
and perturbed rhat"""
def rhat_comps(file_path , file_path_comp = None , pert = None , material = None):
    """input: x1 (array), RK4(5) x vals
              y1 (array), RK4(5) y vals
              t (tuple), consisting of dt and t_tot, time of simulations
              r_per (array), optional, perturbed rhat
              x2 (array), optional, default:None, Leapfrog x vals
              y2 (array), optional, default:None, Leapfrog y vals

        returns: none"""
    
    res = np.load(file_path)
    # x , y , _ , _ , _ , _ , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    x , y , _ , t = [res[k] for k in ("x" , "y" , "b", "t")]
    r = np.sqrt(x**2 + y**2) #r hat
    num_orbits = int(t[-1] / (2 * np.pi))
    t_index = int(len(t) / num_orbits)
    t_numtot = t[::t_index]
    r_numtot = r[::t_index]
    r_betw = np.zeros((len(r_numtot) - 1 , t_index))
    amp = np.zeros((len(r_betw)))

    for i in range(len(r_numtot) - 1):
        r_betw[i, :] = r[i*t_index:(i+1)*t_index]
    
    for i , row in enumerate(r_betw):
        max_ind = np.argmax(i)
        min_ind = np.argmin(i)

        amp[i] = (np.abs(row[max_ind]) + np.abs(row[min_ind])) / 2

    # print(amp)  

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_r.png"

    """comparing RK4(5) and Leapfrog rhat"""
    if file_path_comp is not None:
        res1 = np.load(file_path_comp)
        x1 , y1 , _ , _ , _ , _ , t = [res1[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]

        base_name1 = os.path.splitext(os.path.basename(file_path_comp))[0]
        save_path1 = f"Plots/{base_name}_r_vs_{base_name1}.png"

        r1 = np.sqrt(x1**2 + y1**2) #r hat
        plt.plot(t[::10] , r[::10] , color = "blue", label = r"RK4(5) $\hat{r}$")
        plt.plot(t[::10] , r1[::10] , color = "red" , linestyle = "--" , label = r"Leapfrog $\hat{r}$")
        plt.title(r"$\hat{r}$ from RK4(5) and Leapfrog solution")
        plt.legend()

        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')


    if pert is not None: #comparing RK4(5) with perturbed rhat
        r_per = pert
        rel_fw_err = np.abs(r - r_per) / np.abs(r)
        # print(rel_fw_err)
        save_path2 = f"Plots/{base_name}_r_vs_perturbed.png"
        # plt.scatter(t_numtot , r_numtot , color = "black" , label = "Orbits")
        plt.plot(t[::10] , r[::10] , color = "blue" , label = r"RK4(5) $\hat{r}$")
        plt.plot(t[::10] , r_per[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{r}$")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{r}$")
        plt.title(r"$\hat{r}$ from RK4(5) silicate and perturbed solution")
        plt.legend()
        plt.savefig(save_path2 , dpi = 300 , bbox_inches = 'tight')

    if material == "silicate":
        plt.plot(t[::10] , r[::10])
        plt.title(r"$\hat{r}$ for silicate, real $\hat{\beta}$")
        # plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
        
    
    if material == "carbon":
        plt.plot(t[::10] , r[::10])
        plt.title(r"$\hat{r}$ for carbon, real $\hat{\beta}$")
        # plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
        
    plt.xlabel(r"$\hat{t}$")
    plt.ylabel(r"$\hat{r}$")
    # plt.ylim(0.05 , 1.1)
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
    x , y , vx , vy , _ , _ , t= [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    r = np.sqrt(x**2 + y**2)
    ymin , ymax = np.where(r >= 0.05) , np.where(r <= 1)

    theta_num = np.atan2(y , x) #thetahat
    theta_num = np.unwrap(theta_num) #avoiding discontinuities

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_vr.png"

    v_r = vx * np.cos(theta_num) + vy * np.sin(theta_num) #cartesian to radial vel
    
    if pert is not None:
        pr = np.load(pert)
        _ , _ , _ , v_per , _ , _ , _ = [pr[k] for k in ("omega" , "r" , "theta" , "vr" , "c0" , "b" , "t")]

        save_path1 = f"Plots/{base_name}_vr_vs_perturbed.png"
        plt.plot(t[::100] , v_r[::100] , color = "blue" , label = r"RK4(5) $\hat{v}$")
        plt.plot(t[::100] , v_per[::100] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{v}$")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_r$")
        # if v_r[ymin][-1] < v_per[ymin][-1]:
        #     plt.ylim(v_r[ymin][-1] , v_r[ymax][0])
        
        # else:
        #     plt.ylim(v_per[ymin][-1] , v_per[ymax][0])

        plt.legend()
        plt.title(r"RK4(5) and perturbed $\hat{v}_r$")
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')
    
    if material == "silicate":
        plt.plot(t[::100] , v_r[::100])
        # plt.ylim(v_r[ymin][-1] , v_r[ymax][0])
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_r$")
        plt.title(r"$\hat{v}_r$ for silicate, real $\hat{\beta}$")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    if material == "carbon":
        plt.plot(t[::100] , v_r[::100])
        # plt.ylim(v_r[ymin][-1] , v_r[ymax][0])
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_r$")
        plt.title(r"$\hat{v}_r$ for carbon, real $\hat{\beta}$")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    plt.show()

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
    ymax , ymin = np.where(r >= 0.05) , np.where(r <= 1)

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_omega.png"
    
    theta_num = np.atan2(y , x) #thetahat RK4(5)
    theta_num = np.unwrap(theta_num) #avoiding discontinuities

    angvel_num = (-vx * np.sin(theta_num) + vy * np.cos(theta_num)) / r    

    if pert is not None:
        angvel = pert
        plt.plot(t[::10] , angvel_num[::10] , color = "blue" , label = "RK4(5)")
        plt.plot(t[::10] , angvel[::10] , color = "red" , linestyle = "--" , label = "Perturbed")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\omega}$")
        plt.title(r"$\hat{\omega}$ RK4(5) vs perturbed solution")

        # if angvel_num[ymax][-1] < angvel[ymax][-1]:
        #     plt.ylim(angvel_num[ymin][0] , angvel[ymax][-1])
        
        # else:
        #     plt.ylim(angvel_num[ymin][0] , angvel[ymax][-1])

        plt.legend()
        save_path = f"Plots/{base_name}_omega_vs_perturbed.png"
        
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    if material == "silicate":
        plt.plot(t[::10] , angvel_num[::10])
        # plt.ylim(angvel_num[ymin][0] , angvel_num[ymax][-1])
        plt.title(r"$\hat{\omega}$ for silicate, real $\hat{\beta}$")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    if material == "carbon":
        plt.plot(t[::10] , angvel_num[::10])
        # plt.ylim(angvel_num[ymin][0] , angvel_num[ymax][-1])
        plt.title(r"$\hat{\omega}$ for carbon, real $\hat{\beta}$")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')

    plt.xlabel(r"$\hat{t}$")
    plt.ylabel(r"$\hat{\omega}$")
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
    _ , _ , _ , _ , _ , b , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path2 = f"Plots/{base_name}_beta.png"

    if fw_err == False and b_per is not None and b_analytical is None:
        b_pert = b_per

        save_path3 = f"Plots/{base_name}_beta_interp_vs_expression.png"
        plt.plot(t[::10] , b_pert[::10] , color = "red" , label = r"$\hat{\beta}$ expression")
        plt.plot(t[::10] , b[::10] , color = "blue" , linestyle = "--" , label = r"Interpolated $\hat{\beta}$ curve")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ from interpolated function and from expression")
        plt.legend()
        plt.savefig(save_path3 , dpi = 300 , bbox_inches = 'tight')
        plt.show()
        

    """comparing betahat from RK4(5) to betahat from perturbed and analytical expression"""
    if fw_err == False and b_per is not None and b_analytical is not None:
        pr = np.load(b_per)
        _ , _ , _ , _ , _ , b_pert , _ = [pr[k] for k in ("omega" , "r" , "theta" , "vr" , "c0" , "b" , "t")]

        save_path = f"Plots/{base_name}_beta_vs_analytical_perturbed.png"
        plt.figure()
        plt.plot(t[::10] , b[::10] , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t[::10] , b_pert[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ from RK4(5) and perturbed solution")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.show()

        plt.figure()
        plt.plot(t[::10] , b[::10] , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t[::10] , b_analytical[::10] , color = "orange" , linestyle = "--" , label = r"Analytical $\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ from RK4(5) and analytical solution")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.show()
        
        plt.figure()
        plt.plot(t[::10] , b[::10] , color = "blue" , label = r"RK4(5) $\hat{\beta}$")
        plt.plot(t[::10] , b_pert[::10] , color = "red" , linestyle = "--" , label = r"Perturbed $\hat{\beta}$")
        plt.plot(t[::10] , b_analytical[::10] , color = "orange" , linestyle = "--" , label = r"Analytical $\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ from RK4(5), perturbed and analytical solution")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.legend()
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
        plt.show()
    
    """comparing relative forward errors in beta hat from RK45-perturbed expression and 
    RK45-analytical expression"""
    if fw_err == True and b_per is not None and b_analytical is not None:
        pr = np.load(b_per)
        _ , _ , _ , _ , _ , b_pert , _ = [pr[k] for k in ("omega" , "r" , "theta" , "vr" , "c0" , "b" , "t")]

        save_path1 = f"Plots/{base_name}_beta_rel_fw_err.png"
        fw_err_RK45_per = np.abs(b - b_pert) / np.abs(b)
        fw_err_RK45_analytical = np.abs(b - b_analytical) / np.abs(b)
        plt.plot(t[::10] , fw_err_RK45_per[::10] , color = "blue" , label = "Rel fw error RK4(5) vs perturbed")
        plt.plot(t[::10] , fw_err_RK45_analytical[::10] , color = "red" , linestyle = "--" , label = "Rel fw error RK4(5) vs analytical")
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"Relative forward error")
        plt.title("Relative forward error, RK4(5) vs perturbed and RK4(5) vs analytical")

        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')
        plt.show()
    
    if material == "silicate":
        plt.plot(t[::10] , b[::10])
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ silicate")
        plt.savefig(save_path2 , dpi = 300 , bbox_inches = 'tight')
        plt.show()
    
    if material == "carbon": 
        plt.plot(t[::10] , b[::10])
        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{\beta}$")
        plt.title(r"$\hat{\beta}$ carbon")
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
    x2 , y2 , vx2 , vy2 , m2 , b_vals2 , _ = [res2[k] for k in ("x","y","vx","vy","m","b", "t")] #unpacking solver2

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
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
        plt.show()
    
    if fw_err == True: #plots RK4(5) sols and relative forward errors between RK4(5) and Leapfrog
        save_path1 = f"Plots/{base_name}_energy_rel_fw_err.png"
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
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')
        plt.show()
        
"""plots general beta curves for silicate and carbon"""
def beta_curves(interp = False , comp = False , sample_pts = False):
    """input: None
    
       returns: None"""
    
    """compare interpolated function with true curve for silicate values"""
    if interp == True: 
        interp = pchip(sil_size , sil_betaval)

        plt.plot(sil_size * 10**(-6) , sil_betaval , linestyle = "--" , label = "True curve")
        plt.plot(sil_size * 10**(-6) , interp(sil_size) , linestyle = "-" , label = "Interpolated function")
        plt.savefig(f"Plots/beta_interpolation_curve_silicate.png" , dpi = 300 , bbox_inches = 'tight')
    
    """comparing real silicate and carbon beta curve"""
    if comp == True:
        plt.xscale("log")
        plt.yscale("log")
        plt.plot(sil_size , sil_betaval , color = "red" , linestyle = "-" , label = "Silicate")
        plt.plot(car_size , car_betaval , color = "blue" , linestyle = "--" , label = "Carbon")
        plt.savefig(f"Plots/beta_interpolation_curve_silicate_carbon.png" , dpi = 300 , bbox_inches = 'tight')
    
    """curve with sampling points"""
    if sample_pts == True:
        plt.xscale("log")
        plt.yscale("log")
        plt.scatter([1.54079 * 10**(-6) , 0.17508 * 10**(-6) , 0.04259 * 10**(-6)] , [0.1235 , 0.8560 , 0.2098] , c = "g")
        plt.scatter([1.54079 * 10**(-6) , 0.17508 * 10**(-6) , 0.04259 * 10**(-6)] , [0.2646 , 3.0589 , 1.6179] , c = "orange")
        plt.plot(sil_size * 10**(-6) , sil_betaval , color = "red" , linestyle = "-" , label = "Silicate")
        plt.plot(car_size * 10**(-6) , car_betaval , color = "blue" , linestyle = "--" , label = "Carbon")
        plt.savefig(f"Plots/beta_curve_sampling_points.png" , dpi = 300 , bbox_inches = 'tight')

    plt.title(r"$\beta$ versus particle size")
    plt.xlabel(r"Particle size (m)")
    plt.ylabel(r"$\beta$")
    plt.legend()
    plt.show()

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

    fig , ax = plt.subplots()

    for mat in material: 
        size = sil_size if mat == "silicate" else car_size
        pr = sil_PR if mat == "silicate" else car_PR

        for sw, vals in tsp_vals[mat].items():
            ax.plot(size , vals,
                     color = cl[sw] , linestyle = ls[mat] , label = f"{sw}")

        label = "PR lifetime"
        ax.plot(size , pr , color = "black" , linestyle = ls[mat] , label = label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(0.1 , 10**8)

        ax.set_xlabel(r"Particle size (m)")
        ax.set_ylabel(r"Lifetime (years)")

        ax.set_title("Combined: PR and sputtering lifetimes")
        ax.legend()
        # fig.savefig(f"Plots/{mat}_PR_sputtering_lifetime_tot.png", dpi = 300 , bbox_inches = 'tight')

    plt.show()
        
def PR_spu_lifetime_separate(limits = False):
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
        pr = sil_PR if mat == "silicate" else car_PR

        for sw, vals in tsp_vals[mat].items():
            ax.plot(size , vals,
                     color = cl[sw] , linestyle = ls[mat] , label = f"{sw}")
            
            if mat == "carbon":
                ax.axvspan(0.01516 * 10**(-6) , 0.54840 * 10**(-6) , color = "blue" , alpha = 0.1)

        label = f"{mat} PR" 
        ax.plot(size , pr , color = "black" , linestyle = ls[mat] , label = label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(0.1 , 10**8)

        ax.set_xlabel(r"Particle size (m)")
        ax.set_ylabel(r"Lifetime (years)")

       
        ax.set_title(f"{mat.capitalize()}: PR and sputtering lifetimes")
        ax.legend()
        fig.savefig(f"Plots/{mat}_PR_sputtering_lifetime_separate.png", dpi = 300 , bbox_inches = 'tight')

    plt.show()
    

def v_theta(file_path , pert = None , material = None):

    res = np.load(file_path)
    x , y , vx , vy , _ , _ , t = [res[k] for k in ("x","y","vx","vy","m","b" , "t")] #unpacking file_path
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    save_path = f"Plots/{base_name}_vtheta.png"

    theta = np.atan2(y , x)
    theta = np.unwrap(theta)

    vtheta = -vx * np.sin(theta) + vy * np.cos(theta)

    r = np.sqrt(x**2 + y**2)
    ymax , ymin = np.where(r >= 0.05) , np.where(r <= 1)

    if pert is not None:
        vtheta_pert = pert
        save_path1 = f"Plots/{base_name}_vtheta_pert_comps.png"

        plt.plot(t[::10] , vtheta_pert[::10] , color = "red" , linestyle = "-" , label = "Perturbed")
        plt.plot(t[::10] , vtheta[::10] , color = "blue" , linestyle = "--" , label = "RK4(5)")

        # if vtheta[ymax][-1] > vtheta_pert[ymax][-1]:
        #     plt.ylim(vtheta[ymin][0] , vtheta[ymax][-1])
        
        # else:
        #     plt.ylim(vtheta_pert[ymin][0] , vtheta_pert[ymax][-1])

        plt.xlabel(r"$\hat{t}$")
        plt.ylabel(r"$\hat{v}_{\theta}$")
        plt.title(r"Perturbed and numerical $\hat{v}_{\theta}$")
        plt.legend()
        plt.savefig(save_path1 , dpi = 300 , bbox_inches = 'tight')
    
    if material == "silicate":
        plt.plot(t[::10] , vtheta[::10])
        plt.ylim(vtheta[ymin][0] , vtheta[ymax][-1])
        plt.title(r"$\hat{v}_{\theta}$ silicate")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
    
    if material == "carbon":
        plt.plot(t[::10] , vtheta[::10])
        plt.ylim(vtheta[ymin][0] , vtheta[ymax][-1])
        plt.title(r"$\hat{v}_{\theta}$ carbon")
        plt.savefig(save_path , dpi = 300 , bbox_inches = 'tight')
    
    plt.xlabel(r"$\hat{t}$")
    plt.ylabel(r"$\hat{v}_{\theta}$")
    plt.show()

if __name__ == "__main__":
    par = dust_properties("silicate" , "slow" , "large")
    file_path = "Files/rk45_t6_large_silicate_slowsw.npz"
    res = np.load(file_path)
    # x , y , vx , vy , m , b , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b" , "t")]
    x , y , b , t = [res[k] for k in ("x" , "y" , "b" , "t")]
    
    p = perturbed_functions(par , t , b , find_k = False)
    c0 = p.C0(p.K)
    om , _ , _ = p.omega(p.K)
    r , _ , _ = p.rad(p.K)
    thetaval = p.theta(p.K)
    betas = p.barr
    vthetapert = om * r
    
    rhat = rhat_comps(file_path , file_path_comp = None , pert = None , material = "silicate")
    # thetahat = thetahat_comps(file_path , file_path_comp = None , pert = thetaval)
    # omegahat = omegahat_comps(file_path , pert = om)
    # betahats = b_plot(file_path , betas)
    # vtheta = v_theta(file_path , vthetapert)

    # beta_curves(interp = False , comp = True)
    # PR_spu_lifetime_separate()

    


    
    
    
    


    

    
    
    
    
    