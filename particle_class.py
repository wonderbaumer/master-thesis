import numpy as np
from config import mhat0 , init_cart_scaled , t5 , t6 , t7 , eps
from leapfrog import leapfrog_algorithm
from scipy_solver import particle_motion , pos_vel , arr_variables
from forces_scaled import tot_acc, sputtering

"""class solving scaled equations of motion for a particle using user-specified numerical solver"""
class particle():
    """attributes:                
       sim_time (tuple), consisting of dt and t_tot
        
       solver (string), user-specified solver for equation of motion. 
                        Default:LEAPFROG, else RK45, RK23 or DOP853

       massloss (boolean), optional. default:True if massloss considered, False if not
                   
       methods: 
       
       pos_vel_calcs(), numerical solution for particle parameters based on input of solver type and 
                        if massloss=True or =False"""
    
    def __init__(self , sim_time , epsilon , solver = "LEAPFROG" , massloss = True):
        self.solver = solver
        self.sim_time = sim_time
        self.massloss = massloss
        self.epsilon = epsilon
    
    """calculates position, velocity and other parameters using different solvers"""
    def pos_vel_calcs(self):
        y0 = np.append(init_cart_scaled, mhat0) #initial values for scipy ivp solver
        
        dt , t_tot = self.sim_time #dt and t_tot unpacking
        
        t_span = (0 , t_tot) #time for simulations
        t_eval = np.arange(0 , t_tot , dt) #setting number of timesteps scipy solver

        state = [0 , t_tot / 1000]
        
        if self.solver == "LEAPFROG" and self.massloss == True:
            pos_and_vel1 = leapfrog_algorithm(init_cart_scaled , tot_acc
                     , self.sim_time , self.epsilon , sputtering) #leapfroging using initial cond
            
        elif self.solver == "LEAPFROG" and self.massloss == False:
            pos_and_vel1 = leapfrog_algorithm(init_cart_scaled , tot_acc
                     , self.sim_time) #leapfroging using initial cond
            
        elif self.solver in ["RK45" , "RK23" , "DOP853"] and self.massloss == True: 
            pos_and_vel = particle_motion(pos_vel , t_span , 
                                          y0 , self.solver , t_eval , state , self.epsilon ,  massloss = True) #specified scipy solver
            
            pos_and_vel1 = arr_variables(pos_and_vel) #array of variables
        
        elif self.solver in ["RK45" , "RK23" , "DOP853"] and self.massloss == False: 
            pos_and_vel = particle_motion(pos_vel , t_span , y0 , self.solver , t_eval , state , self.epsilon, massloss = False) #specified scipy solver
            
            pos_and_vel1 = arr_variables(pos_and_vel) #array of variables

        return pos_and_vel1

if __name__ == "__main__":
    
    epsilon = eps(sw = "slow" , species = "all")
    p = particle(t5 , epsilon , "LEAPFROG" , massloss = True)
    vals = p.pos_vel_calcs()
    x , y , vx , vy , m , b = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5]

    #np.savez("Files/rk45_t5_masslossTrue_scaledeqs_dt1.npz" , x = x , y = y , vx = vx , vy = vy , m = m , b = b)