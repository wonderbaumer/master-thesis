from config import *
import sys
sys.path.insert(1, 'C:/Users/Cecilie.Bamer/Documents/Project-paper/')
from particle_class import *
from plot import *
from analytical_functions import *

"""point of execution for simulations and plotting"""
def main(comp_type , time , solver = None , filename1 = None , filename2 = None):
    """input: comp_type (string), eps_beta, theta, r , betahat , energy
              time (tuple), consiting of dt, t_tot, t4 for 500 orbits
              solver (string), default: None, can be LEAPFROG or RK45
              filename1 (string), default: None, name of npz file to open and unpack, numerical solver results
              filename2 (string), default: None, other numerical solver results file"""

    """if filename1 is given, open the file, else run the particle class to produce a file"""
    if filename1 is not None:
        f1 = np.load(filename1)
        solver1 = (f1["x"] , f1["y"] , f1["vx"] , f1["vy"] , f1["m"] , f1["b"])
    else:
        p1 = particle(init_polar , time , solver)
        p1.save_to_file()
        filename1 = f"C:/Users/Cecilie.Bamer/Documents/Project-paper/Files/{p1.solver}_{p1.sim_label}_massloss{p1.massloss}.npz"
        f1 = np.load(filename1)
        solver1 = (f1["x"], f1["y"], f1["vx"], f1["vy"], f1["m"], f1["b"])

    """optional to provide file of solver 2 to compare with solver 1"""   
    if filename2 is not None:
        f2 = np.load(filename2)
        solver2 = (f2["x"] , f2["y"] , f2["vx"] , f2["vy"] , f2["m"] , f2["b"])

    dt , t_tot = time #time unpacking
    t_hat = np.arange(0 , t_tot , dt) / T #hatted time    

    """running plotting codes based on type of comparison user want to make"""
    if comp_type == "eps_beta":
        eps_init_beta() #for epsilon vs B plot for range of masses
    
    elif comp_type == "theta":
        if filename2 is not None: #compare theta for two numerical solvers
            ang_comps(solver1 , time , solver2)

        else: #compare theta with a numerical and perturbed solution
            ang_comps(solver1 , time , angular_position(t_hat))
    
    elif comp_type == "r":
        if filename2 is not None: #compare r two numerical solvers
            rad_comps(solver1 , time , solver2)
        
        else: #compare r one numerical solver and one perturbed solution
            rad_comps(solver1 , time , radial_position(t_hat))
    
    elif comp_type == "betahat":
        #Choose to change fw_err to True to see the relative forward errors num-pert or num-analytical
        b_plot(solver1 , betahat_pert(t_hat) , betahat_analytical(t_hat) , time , fw_err = False)

    elif comp_type == "energy": #compare energy between two numerical solvers
        energy_plot(time , solver1 , solver2)

if __name__ == "__main__":
    """Example of running code"""
    comp_type = "betahat"
    time = t4
    """file path needs to be specified by user"""
    filename1 = "C:/Users/Cecilie.Bamer/Documents/Files/rk45_500orbits_massloss.npz" 

    main(comp_type , time , filename1 = filename1)