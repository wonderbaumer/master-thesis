from config import *
from particle_class import particle
from plot import eps_init_beta , thetahat_comps , rhat_comps , vhat_comps , omegahat_comps , b_plot , energy_plot
from pert_functions import thetahat_pert , rhat_pert , betahat_pert , betahat_analytical , vrhat_pert , omegahat_pert , perturbed_orbit

"""point of execution for simulations and plotting. user can specify comparison type for plotting,
solver type or already saved files to open, can open both files to compare or just one"""
def main(comp_type , time , solver = None , massloss = True , rk_file = None , lf_file = None 
         , fw_err = False):
    """input: comp_type (string), eps_beta, theta, r , betahat , energy
              time (tuple), consiting of dt, t_tot, t4 for 500 orbits
              solver (string), default: None, can be LEAPFROG or RK4(5)
              filename1 (string), default: None, name of npz file to open and unpack, numerical solver results
              filename2 (string), default: None, other numerical solver results file"""
    
    if solver is not None and massloss == True:
        p = particle(sim_time = time, solver = solver , massloss = True) #initializing particle class
        vals = p.pos_vel_calcs() #running numerical solver
        x , y , vx , vy , m , b = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5]
    
    elif solver is not None and massloss == False:
        p = particle(sim_time = time, solver = solver , massloss = False) #initializing particle class
        vals = p.pos_vel_calcs() #running numerical solver
        x , y , vx , vy , m , b = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5]

    """if filename1 is given, open the file, else run the particle class to produce a file"""
    if rk_file is not None:
        f1 = np.load(rk_file)
        x1 , y1 , vx1 , vy1 , m1 , b1 = [f1[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
        rk_vals = (x1 , y1 , vx1 , vy1 , m1 , b1)

    """optional to provide file of solver 2 to compare with solver 1"""   
    if lf_file is not None:
        f2 = np.load(lf_file)
        x2 , y2 , vx2 , vy2 , m2 , b2 = [f2[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]
        lf_vals = (x2 , y2 , vx2 , vy2 , m2, b2)

    dt , t_tot = time #time unpacking
    t_hat = np.arange(0 , t_tot , dt)  #that

    """running plotting codes based on type of comparison user want to make"""
    if comp_type == "eps_beta":
        eps_init_beta() #for epsilon vs B plot for range of masses, only need mandatory inputs to main
    
    elif comp_type == "thetahat": #input rkfile and lffile, only one of the files or a specified solver
        if rk_file is not None and lf_file is not None: #compare thetahat RK4(5) and Leapfrog
            thetahat_comps(x1 , y1 , time , x2 , y2)

        elif rk_file is not None: #compare thetahat RK4(5) and perturbed solution
            thetahat_comps(x1 , y1 , time , theta_per = thetahat_pert(t_hat))
        
        elif lf_file is not None: #compare thetahat Leapfrog and perturbed solution
            thetahat_comps(x2 , y2 , time , theta_per = thetahat_pert(t_hat))
        
        else: #compare thetahat numerical solver and perturbed solution
            thetahat_comps(x , y , time , theta_per = thetahat_pert(t_hat))
    
    elif comp_type == "rhat": #input rkfile and lffile, only one of the files or a specified solver
        if rk_file is not None and lf_file is not None: #compare rhat RK4(5) and Leapfrog
            rhat_comps(x1 , y1 , time , x2 = x2 , y2 = y2)
        
        elif rk_file is not None: #compare rhat RK4(5) and perturbed solution
            rhat_comps(x1 , y1 , time , r_per = rhat_pert(t_hat))
        
        elif lf_file is not None: #compare rhat Leapfrog and perturbed solution
            rhat_comps(x2 , y2 , time , r_per = rhat_pert(t_hat))
        
        else: #compare rhat RK4(5) and one perturbed solution
            rhat_comps(x , y , time , r_per = rhat_pert(t_hat))
    
    elif comp_type == "betahat": #input rkfile, lffile or specified solver, fw_err True or False 
        if rk_file is not None and fw_err == False: #compare betahat RK4(5)
            b_plot(b1 , time , betahat_pert(t_hat) , betahat_analytical(t_hat) , fw_err = False)
        elif rk_file is not None and fw_err == True: #fw err betahat RK4(5)
            b_plot(b1 , time , betahat_pert(t_hat) , betahat_analytical(t_hat) , fw_err = True)
        elif lf_file is not None and fw_err == False: #compare betahat Leapfrog
            b_plot(b2 , time , betahat_pert(t_hat) , betahat_analytical(t_hat) , fw_err = False)
        elif lf_file is not None and fw_err == True: #fw err betahat Leapfrog
            b_plot(b2 , time , betahat_pert(t_hat) , betahat_analytical(t_hat) , fw_err = True)
        elif solver is not None and fw_err == False: #compare betahat numerical solver
            b_plot(b , time , betahat_pert(t_hat) , betahat_analytical(t_hat) , fw_err = False)
        elif solver is not None and fw_err == True: #compare betahat numerical solver
            b_plot(b , time , betahat_pert(t_hat) , betahat_analytical(t_hat) , fw_err = True)


    elif comp_type == "energy": #input rkfile and lffile or one file and specified solver, err True or False

        if rk_file is not None and lf_file is not None and fw_err == False: #compare energy RK4(5) and Leapfrog
            energy_plot(rk_vals , time , lf_vals)
        elif rk_file is not None and lf_file is not None and fw_err == True: #fw err energy RK4(5) and Leapfrog
            energy_plot(rk_vals , time , lf_vals , fw_err = True)
        elif rk_file is not None and solver is not None and fw_err == False: #compare energy RK4(5) and specified solver
            solver_vals = (x , y , vx , vy , m , b)
            energy_plot(rk_vals , time , solver_vals)
        elif rk_file is not None and solver is not None and fw_err == True: #fw err energy RK4(5) and specified solver
            solver_vals = (x , y , vx , vy , m , b)
            energy_plot(rk_vals , time , solver_vals , fw_err = True)
        elif lf_file is not None and solver is not None and fw_err == False: #compare energy Leapfrog and specified solver
            solver_vals = (x , y , vx , vy , m , b)
            energy_plot(lf_vals , time , solver_vals)
        elif lf_file is not None and solver is not None and fw_err == True: #fw err energy Leapfrog and specified solver
            solver_vals = (x , y , vx , vy , m , b)
            energy_plot(lf_vals , time , solver_vals , fw_err = True)
    
    elif comp_type == "vhat": #input rkfile, lffile or specified solver
        v_pert = vrhat_pert(t_hat)

        if rk_file is not None:
            vhat_comps(x1 , y1 , vx1 , vy1 , time , v_pert)
        
        elif lf_file is not None:
            vhat_comps(x2 , y2 , vx2 , vy2 , time , v_pert)
        
        elif solver is not None:
            vhat_comps(x , y , vx , vy , time , v_pert)
    
    elif comp_type == "omegahat":
        ang_vel = omegahat_pert(t_hat)

        if rk_file is not None:
            omegahat_comps(x1 , x1 , vx1 , vy1 , time , ang_vel)
        
        elif lf_file is not None:
            omegahat_comps(x2 , y2 , vx2 , vy2 , time , ang_vel)
        
        elif solver is not None:
            omegahat_comps(x , y , vx , vy , time , ang_vel)


if __name__ == "__main__":
    """Example of running code"""
    comp_type = "vhat"
    rk_file = "Files/rk45_t6_masslossTrue_scaledeqs.npz"
    lf_file = "Files/leapfrog_t6_masslossTrue_scaledeqs.npz"
    

    main(comp_type , t6 , solver = "RK45" , rk_file = None , lf_file = None , massloss = True , fw_err = False)