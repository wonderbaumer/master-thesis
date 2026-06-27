import numpy as np
from config import t6 , t7 , t8 , t9 , t10
from particle_class import particle_solver
from plot import (eps_init_betareal , rhat_comps , omegahat_comps , b_plot 
                  , energy_plot , beta_curves , PR_spu_lifetime_num , eval_sizes , ecc_sc , ecc_math 
                  , PR_spu_lifetime_pert , PR_spu_lifetime_theo)
from pert_variable_eps import runner_class
from new_num_lifetimes import true_lifetime , true_lifetime_variableeps
from dust_properties import dust_properties

"""Function plotting all plots displayed in master's thesis"""
def main(plot_type , particle_obj = None , pert = None , material = "silicate" , lf_path = None 
         , rk_path = None , mat_comp = None , zoomed = False):
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
              mat_comp (string), default:None, else path to comparison .npz file
              zoomed (bool), default:False, else given plot scaled to preselected values
              """
    
    if plot_type == "exp_beta":
        beta_curves(pert = True)
    
    if plot_type == "expected_lifetimes":
        PR_spu_lifetime_theo()
    
    if plot_type == "eps_vs_B":
        eps_init_betareal()
    
    if plot_type == "beta_interp":
        beta_curves(True , "silicate")
    
    if plot_type == "init_sizes":
        eval_sizes()
    
    if plot_type == "energy_comp":
        energy_plot(rk_path , lf_path , particle_obj , False)
    
    if plot_type == "beta_comp":
        b_plot(rk_path , pert , material , mat_comp , zoomed)

    if plot_type == "r_comp":
        rhat_comps(rk_path , pert , material , mat_comp , zoomed)
    
    if plot_type == "omega_comp":
        omegahat_comps(rk_path , pert , material , mat_comp)
    
    if plot_type == "lifetimes_comp":
        PR_spu_lifetime_pert()
    
    if plot_type == "ecc_scaled":
        ecc_sc(rk_path , particle_obj , None)
    
    if plot_type == "ecc_math":
        ecc_math(rk_path , None)

def make_file(mat , swcond , parsize , sim_time , filename , massloss , drag , solver , analytical = False):
    par = dust_properties(material = mat , sw = swcond , size = parsize)
    p = particle_solver(sim_time , par , solver , massloss , drag , analytical)
    vals = p.pos_vel_calcs()
    x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6] 
    np.savez(f"Files/{filename}.npz" , x = x[::10] , y = y[::10] , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10]) 

"""Neccessary files"""
# file1 = make_file("silicate" , "slow" , "A" , t6 , "rk45_t6_A_silicate_slowsw_nomassloss_nodrag" 
# , False , False , "RK45")

# file2 = make_file("silicate" , "slow" , "A" , t6 , "leapfrog_t6_A_silicate_slowsw_nomassloss_nodrag" 
# , False , False , "LEAPFROG")

# file3 = make_file("silicate" , "slow" , "C" , t8 , "rk45_t8_C_silicate_slowsw" , True , True 
#                   , "RK45" , True)

# file4 = make_file("carbon" , "slow" , "C" , t8 , "rk45_t8_C_carbon_slowsw" , True , True 
# , "RK45" , True)

# file5 = make_file("carbon" , "slow" , "C" , t8 , "rk45_t8_C_carbon_slowsw_numbeta" , True 
# , True , "RK45")

# file6 = make_file("silicate" , "slow" , "C" , t8 , "rk45_t8_C_silicate_slowsw_numbeta" , True 
# , True , "RK45")
         
# file7 = make_file("silicate" , "slow" , "A" , t8 , "rk45_t8_A_silicate_slowsw_numbeta" , True 
# , True , "RK45" , False)

# file8 = make_file("silicate" , "CME" , "C" , t8 , "rk45_t8_C_silicate_CMEsw_numbeta" , True 
# , True , "RK45")

# file9 = make_file("silicate" , "slow" , "C" , t6 , "rk45_t6_C_silicate_slowsw_numbeta" , True 
# , True , "RK45" , False)

"""numerical and analytical beta curve"""
# main("exp_beta")

"""theoretical lifetimes"""
# main("expected_lifetimes")

"""epsilon_vs_B"""
# main("eps_vs_B")

"""interpolated beta curve"""
# main("beta_interp")

"""sampling points"""
# main("init_sizes")

"""energy comparisons , file 1 and 2"""
# par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "A")

# main("energy_comp" , particle_obj = par 
#      , rk_path = "Files/rk45_t6_A_silicate_slowsw_nomassloss_nodrag.npz" 
#      , lf_path = "Files/leapfrog_t6_A_silicate_slowsw_nomassloss_nodrag.npz")

"""multiscale validation, silicate r, all analytical beta , file 3"""
# par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "C")
         
# p1 = runner_class(par , t8)
# vals1 = p1.solver() 
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# rpert = (r0 , tpert / par.epsilon)

# main("r_comp" , rk_path = "Files/rk45_t8_C_silicate_slowsw.npz" , pert = rpert , material = "silicate")

"""multiscale validation, silicate omega, all analytical beta , file 3"""  
# par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "C")
#        
# p1 = runner_class(par , t8)
# vals1 = p1.solver()  
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# omegapert = (omega0 , tpert / par.epsilon)

# main("omega_comp" , rk_path = "Files/rk45_t8_C_silicate_slowsw.npz" , pert = omegapert 
#      , material = "silicate")

"""multiscale validation, silicate beta, all analytical beta, file 3"""
# par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "C")
         
# p1 = runner_class(par , t8)
# vals1 = p1.solver()
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# bpert = (beta0 , tpert / par.epsilon)

# main("beta_comp" , rk_path = "Files/rk45_t8_C_silicate_slowsw.npz" , pert = bpert 
#      , material = "silicate")

"""multiscale validation, carbon r, all analytical beta , file 4"""
# par = dust_properties("carbon" , "slow" , init_dist = 1.0 , size = "C")
         
# p1 = runner_class(par , t8)
# vals1 = p1.solver() 
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# rpert = (r0 , tpert / par.epsilon)

# main("r_comp" , rk_path = "Files/rk45_t8_C_carbon_slowsw.npz" , pert = rpert , material = "carbon")

"""multiscale validation, carbon omega, all analytical beta, file 4"""
# par = dust_properties("carbon" , "slow" , init_dist = 1.0 , size = "C")

# p1 = runner_class(par , t8)
# vals1 = p1.solver() 
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# omegapert = (omega0 , tpert / par.epsilon)

# main("omega_comp" , rk_path = "Files/rk45_t8_C_carbon_slowsw.npz" , pert = omegapert 
#      , material = "carbon")

"""multiscale validation, carbon beta, all analytical beta, file 4"""
# par = dust_properties("carbon" , "slow" , init_dist = 1.0 , size = "C")

# p1 = runner_class(par , t8)
# vals1 = p1.solver()   
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# bpert = (beta0 , tpert / par.epsilon)

# main("beta_comp" , rk_path = "Files/rk45_t8_C_carbon_slowsw.npz" , pert = bpert 
#      , material = "carbon")

"""numerical beta validation, silicate and carbon r, file 5 and 6"""
# main("r_comp" , rk_path = "Files/rk45_t8_C_silicate_slowsw_numbeta.npz" , mat_comp = "Files/rk45_t8_C_carbon_slowsw_numbeta.npz")

"""numerical beta validation, silicate and carbon omega, file 5 and 6"""
# main("omega_comp" , rk_path = "Files/rk45_t8_C_silicate_slowsw_numbeta.npz" , mat_comp = "Files/rk45_t8_C_carbon_slowsw_numbeta.npz")

"""numerical beta validation, silicate and carbon beta, file 5 and 6"""
# main("beta_comp" , rk_path = "Files/rk45_t8_C_silicate_slowsw_numbeta.npz" , mat_comp = "Files/rk45_t8_C_carbon_slowsw_numbeta.npz")

"""numerical beta validation, silicate size A, r, file 7"""
# main("r_comp" , rk_path = "Files/rk45_t8_A_silicate_slowsw_numbeta.npz" , material = "silicate")

"""numerical beta validation, silicate size A, omega, file 7"""
# main("omega_comp" , rk_path = "Files/rk45_t8_A_silicate_slowsw_numbeta.npz" , material = "silicate")


"""numerical beta validation, silicate size A, beta, file 7"""      
# main("beta_comp" , rk_path = "Files/rk45_t8_A_silicate_slowsw_numbeta.npz" , material = "silicate")

"""comparison analytical and numerical beta, silicate r, file 6"""
# par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "C")

# p1 = runner_class(par , t8)
# vals1 = p1.solver()
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# rpert = (r0 , tpert / par.epsilon)

# main("r_comp" , rk_path = "Files/rk45_t8_C_silicate_slowsw_numbeta.npz" , pert = rpert , material = "silicate")

"""comparison analytical and numerical beta, silicate omega, file 6"""
# par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "C")

# p1 = runner_class(par , t8)
# vals1 = p1.solver()
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# omegapert = (omega0 , tpert / par.epsilon)

# main("omega_comp" , rk_path = "Files/rk45_t8_C_silicate_slowsw_numbeta.npz" , pert = omegapert , material = "silicate")

"""comparison analytical and numerical beta, silicate beta, file 6"""
# par = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "C")

# p1 = runner_class(par , t8)
# vals1 = p1.solver()
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# bpert = (beta0 , tpert / par.epsilon)

# main("beta_comp" , rk_path = "Files/rk45_t8_C_silicate_slowsw_numbeta.npz" , pert = bpert , material = "silicate")

"""comparison analytical and numerical beta, carbon r, file 5"""
# par = dust_properties("carbon" , "slow" , init_dist = 1.0 , size = "C")

# p1 = runner_class(par , t8)
# vals1 = p1.solver()
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# rpert = (r0 , tpert / par.epsilon)

# main("r_comp" , rk_path = "Files/rk45_t8_C_carbon_slowsw_numbeta.npz" , pert = rpert , material = "carbon")

"""comparison analytical and numerical beta, carbon omega, file 5"""
# par = dust_properties("carbon" , "slow" , init_dist = 1.0 , size = "C")

# p1 = runner_class(par , t8)
# vals1 = p1.solver()
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# omegapert = (omega0 , tpert / par.epsilon)

# main("omega_comp" , rk_path = "Files/rk45_t8_C_carbon_slowsw_numbeta.npz" , pert = omegapert , material = "carbon")

"""comparison analytical and numerical beta, carbon beta, file 5"""
# par = dust_properties("carbon" , "slow" , init_dist = 1.0 , size = "C")

# p1 = runner_class(par , t8)
# vals1 = p1.solver()
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# bpert = (beta0 , tpert / par.epsilon)

# main("beta_comp" , rk_path = "Files/rk45_t8_C_carbon_slowsw_numbeta.npz" , pert = bpert , material = "carbon")

"""numerical and perturbed vs theoretical lifetimes"""
# main("lifetimes_comp")

"""numerical and perturbed silicate particle C in CME conditions, r, file 8"""
# parsil = dust_properties("silicate" , "CME" , init_dist = 1.0 , size = "C")

# p1 = runner_class(parsil , t8)
# vals1 = p1.solver()
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# rpert = (r0 , tpert / parsil.epsilon)
         
# main("r_comp" , rk_path = "Files/rk45_t8_C_silicate_CMEsw_numbeta.npz"
#      , material = "silicate" , zoomed = True , pert = rpert)

"""numerical and perturbed silicate particle C in CME conditions, beta, file 8"""
# parsil = dust_properties("silicate" , "CME" , init_dist = 1.0 , size = "C")

# p1 = runner_class(parsil , t8)
# vals1 = p1.solver()
# r0 , omega0 , beta0 , tpert , C0 = vals1[0] , vals1[1] , vals1[2] , vals1[3] , vals1[4]
# bpert = (beta0 , tpert / parsil.epsilon)
         
# main("beta_comp" , rk_path = "Files/rk45_t8_C_silicate_CMEsw_numbeta.npz"
#      , material = "silicate" , zoomed = True , pert = bpert)

"""eccentricity calcs, file 9"""
# parsil = dust_properties("silicate" , "slow" , init_dist = 1.0 , size = "C")

# main("ecc_math" , rk_path = "Files/rk45_t6_C_silicate_slowsw_numbeta.npz")
# main("ecc_scaled" , rk_path = "Files/rk45_t6_C_silicate_slowsw_numbeta.npz" , particle_obj = parsil)