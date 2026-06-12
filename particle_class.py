import numpy as np
from leapfrog import leapfrog_algorithm
from scipy_solver import particle_motion , pos_vel , arr_variables
from forces_scaled import tot_acc , sputtering , inter_func
from config import t6 , t7 , t8 , t9 , t10 , material_files , betahat0 , mhat0 , init_cart_scaled
from dust_properties import dust_properties


class particle_solver():
    """Class solving scaled equations of motion for a particle using user-specified numerical solver
    
    Attributes:                
       sim_time (tuple), consisting of dt and t_tot
       par (class), particle properties based on initial distance, radius and B, solar wind conds
                    and material
       solver (string), user-specified solver for equation of motion. 
                        Default:RK45, else LEAPFROG, RK23 or DOP853
       massloss (boolean), optional. default:True if massloss considered, False if not

       r (float), initial distance of particle
       B (float), initial beta value
       V (float), initial orbital velocity
       T (float), initial orbital period
       delta (float), drag term
       epsilon (float), mass loss rate
       K (float), drag-to-mass loss ratio
       mhat0 (float), initial scaled mass
       betahat0 (float), initial scaled beta
       init_cart_scaled (array), initial scaled cartesian coordinates
       file (.dat), r and beta file corresponding to input material
       beta_func (function), interpolated beta(r) from file values
       

                   
    Methods: 
       
       pos_vel_calcs(), numerical solution for particle parameters based on input of solver type and 
                        if massloss=True or False"""
    
    def __init__(self , sim_time , par , solver = "RK45" , massloss = True):
        self.solver = solver
        self.sim_time = sim_time
        self.massloss = massloss
        self.particle = par
        self.r = par.r
        self.B = par.B
        self.V = par.V
        self.T = par.T
        self.delta = par.delta
        self.epsilon = par.epsilon
        self.K = par.K
        self.mhat0 = mhat0
        self.betahat0 = betahat0
        self.init_cart_scaled = init_cart_scaled
        self.file = material_files[self.particle.material]
        self.beta_func = inter_func(self.file) 

    """Calculates position, velocity, mass and beta with user-specified numerical solver"""
    def pos_vel_calcs(self):
        y0 = np.concatenate((self.init_cart_scaled, [self.mhat0])) #initial scaled cartesian values
        
        dt , t_tot = self.sim_time #dt and t_tot unpacking
        dt = dt  
        t_tot = t_tot 
        t_span = (0 , t_tot) #time for simulations
        state = [0 , t_tot / 1000] #for progress bar
        
        stopping_reason = None
        stopping_events = ["particle outside interpolation range" , "particle impacted Sun"]
        
        """Simulations based on solver and mass loss specifications"""
        if self.solver == "LEAPFROG" and self.massloss == True:
            pos_and_vel1 = leapfrog_algorithm(y0 , tot_acc , self.sim_time , self , self.epsilon 
                                              , sputtering) 
            
        elif self.solver == "LEAPFROG" and self.massloss == False:
            pos_and_vel1 = leapfrog_algorithm(y0 , tot_acc , self.sim_time , self)
            
        elif self.solver in ["Radau" , "RK45" , "RK23" , "DOP853"] and self.massloss == True: 
            pos_and_vel = particle_motion(pos_vel , t_span , y0 , self.solver , state , self.epsilon ,
                                          self , massloss = True)
            
            for i , te in enumerate(pos_and_vel.t_events):
                if len(te) > 0:
                    stopping_reason = stopping_events[i]
                    break

            print(stopping_reason)

            pos_and_vel1 = arr_variables(pos_and_vel , self , self.epsilon , True)
        
        elif self.solver in ["RK45" , "RK23" , "DOP853"] and self.massloss == False: 
            pos_and_vel = particle_motion(pos_vel , t_span , y0 , self.solver , state , self.epsilon ,
                                          self , massloss = False)

            for i , te in enumerate(pos_and_vel.t_events):
                if len(te) > 0:
                    stopping_reason = stopping_events[i]
                    break

            print(stopping_reason)

            pos_and_vel1 = arr_variables(pos_and_vel , self , self.epsilon , False)

        return pos_and_vel1

if __name__ == "__main__":
    par = dust_properties("silicate" , "slow" , init_dist = 1 , size = "E")
    p = particle_solver(t6 , par , "RK45" , massloss = True)
    vals = p.pos_vel_calcs()
    
    x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6] , vals[: , 7]
    
    # np.savez("Files/rk45_t6_E_silicate_slowsw.npz" , x = x[::10] , y = y[::10] , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10])
    
    

    
    