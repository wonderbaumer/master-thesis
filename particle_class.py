import numpy as np
from config import mhat0 , init_cart_scaled , t4 , t5 , t6 , t7
from leapfrog import leapfrog_algorithm
from scipy_solver import particle_motion , pos_vel , arr_variables
from forces_scaled import tot_acc, sputtering

"""class solving scaled equations of motion for a particle using user-specified numerical solver"""
class particle():
    """attributes:
       
       init_cond (array), initial scaled x , y , vx , vy
                          
       sim_time (tuple), consisting of dt and t_tot
        
       solver (string), user-specified solver for equation of motion. 
                        Default:LEAPFROG, else RK45, RK23 or DOP853

       massloss (boolean), optional. default:True if massloss considered, False if not
                   
       methods: 
       
       pos_vel_calcs(), numerical solution for particle parameters based on input of solver type and 
                        if massloss=True or =False"""
    
    def __init__(self , init_cond , sim_time , solver = "LEAPFROG" , massloss = True):
        self.init_cond = init_cond
        self.solver = solver
        self.sim_time = sim_time
        self.massloss = massloss
    
    """calculates position, velocity and other parameters using different solvers"""
    def pos_vel_calcs(self):
        initial_vals = self.init_cond #initial values for leapfrog
        y0 = np.append(initial_vals , mhat0) #initial values for scipy ivp solver
        
        dt , t_tot = self.sim_time #dt and t_tot unpacking
        
        t_span = (0 , t_tot) #time for simulations
        t_eval = np.arange(0 , t_tot , dt) #setting number of timesteps scipy solver
        
        if self.solver == "LEAPFROG" and self.massloss == True:
            pos_and_vel1 = leapfrog_algorithm(initial_vals , tot_acc
                     , self.sim_time , sputtering) #leapfroging using initial cond
            
        elif self.solver == "LEAPFROG" and self.massloss == False:
            pos_and_vel1 = leapfrog_algorithm(initial_vals , tot_acc
                     , self.sim_time) #leapfroging using initial cond
            
        elif self.solver in ["RK45" , "RK23" , "DOP853"] and self.massloss == True: 
            pos_and_vel = particle_motion(pos_vel , t_span , 
                                          y0 , self.solver , t_eval , massloss = True) #specified scipy solver
            
            pos_and_vel1 = arr_variables(pos_and_vel) #array of variables
        
        elif self.solver in ["RK45" , "RK23" , "DOP853"] and self.massloss == False: 
            pos_and_vel = particle_motion(pos_vel , t_span , y0 , self.solver , t_eval
                                          , massloss = False) #specified scipy solver
            
            pos_and_vel1 = arr_variables(pos_and_vel) #array of variables

        return pos_and_vel1

if __name__ == "__main__":
    
    p = particle(init_cart_scaled , t6 , "LEAPFROG" , massloss = False)
    vals = p.pos_vel_calcs()
    x , y , vx , vy , m , b = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5]
    
    print(np.sqrt(x**2 + y**2))