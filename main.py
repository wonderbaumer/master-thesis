import numpy as np
from config import t6 , t7 , t8 , t9 , t10
from particle_class import particle_solver
from plot import (eps_init_betareal , rhat_comps , omegahat_comps , b_plot 
                  , energy_plot , beta_curves , PR_spu_lifetime , eval_sizes , ecc_sc , ecc_math 
                  , PR_spu_lifetime_separate)
from pert_variable_eps import rhs , pert_motion
from lifetime_calcs import true_lifetime_variableeps , true_lifetime
from dust_properties import dust_properties

"""Function plotting all plots displayed in master's thesis"""
def main(plot_type , particle_obj = None , pert = None , material = "silicate" , lf_path = None 
         , rk_path = None):
    """input: plot_type (string), options: 
              beta_curves, experimental beta curves for silicate and carbon
              expected_lifetimes, theoretical pr and sputtering lifetimes both materials
              eps_vs_B, epsilon for all B values, both materials
              beta_interp, interpolated beta curve for silicate
              eval_sizes, initial particle sizes on beta curves
              beta_comp, rk4(5) and pert beta comparison
              r_comp, rk4(5) and pert r comparison
              omega_comp, rk4(5) and pert omega comparison
              lifetimes_cst_eps, rk4(5) vs theoretical lifetimes with constant sputtering
              lifetimes_varying_eps, rk4(5) vs theoretical lifetime with varying sputtering
              num_r, rk4(5) r only
              num_beta, rk4(5) beta only
              num_omega, rk4(5) omega only
              ecc_scaled, rk4(5) eccentricity from scaled expression
              ecc_math, rk4(5) eccentricity from elliptical definitions

              particle_obj (instance), containing particle properties
              pert (array), default:None, else relevant perturbed parameter
              material (string), default: silicate, else carbon
              lf_path (string), default:None, else path to leapfrog file
              rk_path (string), default:None, else path to rk4(5) file
              """
    
    if plot_type == "exp_beta":
        beta_curves(comp = True)
    
    if plot_type == "expected_lifetimes":
        PR_spu_lifetime()
    
    if plot_type == "eps_vs_B":
        eps_init_betareal()
    
    if plot_type == "beta_interp":
        beta_curves(True , "silicate")
    
    if plot_type == "init_sizes":
        eval_sizes()
    
    if plot_type == "energy_comp":
        energy_plot(rk_path , lf_path , particle_obj , False)
    
    if plot_type == "beta_comp":
        b_plot(rk_path , pert , material)

    if plot_type == "r_comp":
        rhat_comps(rk_path , pert , material)
    
    if plot_type == "omega_comp":
        omegahat_comps(rk_path , pert , material)
    
    if plot_type == "lifetimes_cst_eps":
        PR_spu_lifetime_separate(file = true_lifetime , lifetime_effects = "both")
    
    if plot_type == "lifetimes_varying_eps":
        PR_spu_lifetime_separate(file = true_lifetime_variableeps , lifetime_effects = "both")
    
    if plot_type == "num_r":
        rhat_comps(rk_path , None , material)
    
    if plot_type == "num_beta":
        b_plot(rk_path , None , material)

    if plot_type == "num_omega":
        omegahat_comps(rk_path , None , material)
    
    if plot_type == "ecc_scaled":
        ecc_sc(rk_path , particle_obj , None)
    
    if plot_type == "ecc_math":
        ecc_math(rk_path , None)

"""experimental beta curve"""
# main("exp_beta")

"""theoretical lifetimes"""
# main("expected_lifetimes")

"""epsilon_vs_B"""
# main("eps_vs_B")

"""interpolated beta curve"""
# main("beta_interp")

"""sampling points"""
# main("init_sizes")

"""energy comparisons"""
# par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "A")
# p = particle_solver(t6 , par , "RK45" , massloss = False , drag = False)
# vals = p.pos_vel_calcs()
    
# x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6] 
    
# np.savez("Files/rk45_t6_A_silicate_slowsw_nomassloss_nodrag.npz" , x = x[::10] , y = y[::10] 
#          , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10]) 

# p = particle_solver(t6 , par , "LEAPFROG" , massloss = False , drag = False)
# vals = p.pos_vel_calcs()
    
# x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6] 
    
# np.savez("Files/leapfrog_t6_A_silicate_slowsw_nomassloss_nodrag.npz" , x = x[::10] , y = y[::10] 
#          , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10]) 

# main("energy_comp" , particle_obj = par 
#      , rk_path = "Files/rk45_t6_A_silicate_slowsw_nomassloss_nodrag.npz" 
#      , lf_path = "Files/leapfrog_t6_A_silicate_slowsw_nomassloss_nodrag.npz")

"""beta curve rk4(5) silicate size A"""
# par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "A")
# p = particle_solver(t6 , par , "RK45" , massloss = True , drag = True)
# vals = p.pos_vel_calcs()
    
# x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6] 
    
# np.savez("Files/rk45_t6_A_silicate_slowsw.npz" , x = x[::10] , y = y[::10] 
#          , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10])

# main("num_beta" , rk_path = "Files/rk45_t6_A_silicate_slowsw.npz" , material = "silicate")

"""beta curve rk4(5) carbon size A"""
# par = dust_properties("carbon" , "slow" , init_dist = 1.0 , size = "A")
# p = particle_solver(t6 , par , "RK45" , massloss = True , drag = True)
# vals = p.pos_vel_calcs()
    
# x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6] 
    
# np.savez("Files/rk45_t6_A_carbon_slowsw.npz" , x = x[::10] , y = y[::10] 
#          , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10])

# main("num_beta" , rk_path = "Files/rk45_t6_A_carbon_slowsw.npz" , material = "carbon")

"""r curve rk4(5) silicate size A"""
# par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "A")
# p = particle_solver(t6 , par , "RK45" , massloss = True , drag = True)
# vals = p.pos_vel_calcs()
    
# x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6] 
    
# np.savez("Files/rk45_t6_A_silicate_slowsw.npz" , x = x[::10] , y = y[::10] 
#          , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10])

# main("num_r" , rk_path = "Files/rk45_t6_A_silicate_slowsw.npz" , material = "silicate")

"""r curve rk4(5) carbon size A"""
# par = dust_properties("carbon" , "slow" , init_dist = 1.0 , size = "A")
# p = particle_solver(t6 , par , "RK45" , massloss = True , drag = True)
# vals = p.pos_vel_calcs()
    
# x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6] 
    
# np.savez("Files/rk45_t6_A_carbon_slowsw.npz" , x = x[::10] , y = y[::10] 
#          , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10])

# main("num_r" , rk_path = "Files/rk45_t6_A_carbon_slowsw.npz" , material = "carbon")

"""omega curve rk4(5) silicate size A"""
# par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "A")
# p = particle_solver(t6 , par , "RK45" , massloss = True , drag = True)
# vals = p.pos_vel_calcs()
    
# x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6] 
    
# np.savez("Files/rk45_t6_A_silicate_slowsw.npz" , x = x[::10] , y = y[::10] 
#          , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10])

# main("num_omega" , rk_path = "Files/rk45_t6_A_silicate_slowsw.npz" , material = "silicate")

"""omega curve rk4(5) carbon size A"""
# par = dust_properties("carbon" , "slow" , init_dist = 1.0 , size = "A")
# p = particle_solver(t6 , par , "RK45" , massloss = True , drag = True)
# vals = p.pos_vel_calcs()
    
# x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6] 
    
# np.savez("Files/rk45_t6_A_carbon_slowsw.npz" , x = x[::10] , y = y[::10] 
#          , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10])

# main("num_omega" , rk_path = "Files/rk45_t6_A_carbon_slowsw.npz" , material = "carbon")

"""beta curve rk4(5) silicate size C"""
# par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "C")
# p = particle_solver(t6 , par , "RK45" , massloss = True , drag = True)
# vals = p.pos_vel_calcs()
    
# x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6] 
    
# np.savez("Files/rk45_t6_C_silicate_slowsw.npz" , x = x[::10] , y = y[::10] 
#          , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10])

# main("num_beta" , rk_path = "Files/rk45_t6_C_silicate_slowsw.npz" , material = "silicate")

"""r curve rk4(5) silicate size C"""
# par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "C")
# p = particle_solver(t6 , par , "RK45" , massloss = True , drag = True)
# vals = p.pos_vel_calcs()
    
# x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6] 
    
# np.savez("Files/rk45_t6_C_silicate_slowsw.npz" , x = x[::10] , y = y[::10] 
#          , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10])

# main("num_r" , rk_path = "Files/rk45_t6_C_silicate_slowsw.npz" , material = "silicate")

"""rk4(5) vs theoretical lifetimes, constant sputtering"""
# main("lifetimes_cst_eps")

"""rk4(5) vs theoretical lifetimes, varying sputtering"""
# main("lifetimes_varying_eps")

"""rk4(5) vs perturbed r, silicate"""
par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "A")
p = particle_solver(t6 , par , "RK45" , massloss = True , drag = True)
vals = p.pos_vel_calcs()
    
x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6] 
    
np.savez("Files/rk45_t6_A_silicate_slowsw.npz" , x = x[::10] , y = y[::10] 
         , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10])

# res = np.load("Files/rk45_t6_A_silicate_slowsw.npz")
# x , y , _ , _ , m , b , t  = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]

func = rhs(t * par.epsilon , [1.0 , 1.0] , par.B , par.K)
r0 , _ , _ , _ = pert_motion(rhs , t[::10] * par.epsilon , [1.0 , 1.0] , par.B , par.K)

main(plot_type = "r_comp" , rk_path = "Files/rk45_t6_A_silicate_slowsw.npz" , pert = r0 , material = "silicate")
