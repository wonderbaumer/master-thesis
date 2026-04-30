import numpy as np
from leapfrog import leapfrog_algorithm
from scipy_solver import particle_motion , pos_vel , arr_variables
from forces_scaled import tot_acc , sputtering , inter_func
from config import t5 , t6 , t7 , material_files , betahat0 , mhat0 , init_cart_scaled , yr
from scipy.constants import G
from polar_to_cart import polar_to_cartesian
from dust_properties import dust_properties
import matplotlib.pyplot as plt

""""""
class particle_solver():
    """class solving scaled equations of motion for a particle using user-specified numerical solver
    
    Attributes:                
       sim_time (tuple), consisting of dt and t_tot
        
       solver (string), user-specified solver for equation of motion. 
                        Default:LEAPFROG, else RK45, RK23 or DOP853

       massloss (boolean), optional. default:True if massloss considered, False if not

                   
       methods: 
       
       pos_vel_calcs(), numerical solution for particle parameters based on input of solver type and 
                        if massloss=True or =False"""
    
    def __init__(self , sim_time , par , solver = "LEAPFROG" , massloss = True):
        self.solver = solver
        self.sim_time = sim_time
        self.massloss = massloss
        self.particle = par
        self.r = par.r
        self.B = par.B
        self.V = par.V
        self.T = par.T
        self.K = par.K
        self.delta = par.delta
        self.epsilon = par.epsilon
        self.mhat0 = mhat0
        self.betahat0 = betahat0

        self.init_cart_scaled = init_cart_scaled
        self.file = material_files[self.particle.material]
        
        self.beta_func = inter_func(self.file) #interpolation function for beta values

    """calculates position, velocity and other parameters using different solvers"""
    def pos_vel_calcs(self):
        y0 = np.concatenate((self.init_cart_scaled, [self.mhat0])) #initial values for scipy ivp solver
        
        dt , t_tot = self.sim_time #dt and t_tot unpacking
        dt = dt / self.T 
        t_tot = t_tot 
        t_span = (0 , t_tot) #time for simulations
        t_eval = np.arange(0 , t_tot , dt) #setting number of timesteps scipy solver
        state = [0 , t_tot / 1000]
        
        stopping_reason = None
        stopping_events = ["particle outside interpolation range" , "particle impacted Sun"]
        
        if self.solver == "LEAPFROG" and self.massloss == True:
            pos_and_vel1 = leapfrog_algorithm(self.init_cart_scaled , tot_acc
                     , self.sim_time , self , self.epsilon , sputtering) #leapfroging using initial cond
            
        elif self.solver == "LEAPFROG" and self.massloss == False:
            pos_and_vel1 = leapfrog_algorithm(self.init_cart_scaled , tot_acc
                     , self.sim_time , self) #leapfroging using initial cond
            
        elif self.solver in ["Radau" , "RK45" , "RK23" , "DOP853"] and self.massloss == True: 
            pos_and_vel = particle_motion(
                pos_vel,
                t_span,
                y0,
                self.solver,
                t_eval,
                state,
                self.epsilon,
                self,              
                massloss=True
                )
            
            for i , te in enumerate(pos_and_vel.t_events):
                if len(te) > 0:
                    stopping_reason = stopping_events[i]
                    break
            print(stopping_reason)

            pos_and_vel1 = arr_variables(pos_and_vel, self)
        
        elif self.solver in ["RK45" , "RK23" , "DOP853"] and self.massloss == False: 
            pos_and_vel = particle_motion(
                pos_vel,
                t_span,
                y0,
                self.solver,
                t_eval,
                state,
                self.epsilon,
                self,              
                massloss=False
            )

            for i , te in enumerate(pos_and_vel.t_events):
                if len(te) > 0:
                    stopping_reason = stopping_events[i]
                    break

            print(stopping_reason)
            pos_and_vel1 = arr_variables(pos_and_vel, self)

        return pos_and_vel1

if __name__ == "__main__":
    par = dust_properties("silicate" , "CME" , "large")
    p = particle_solver(t7 , par , "RK45" , massloss = True)
    vals = p.pos_vel_calcs()

    x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6]
    print(t * p.T / yr , np.sqrt(x**2 + y**2)[-1] , m)
    # np.savez("Files/rk45_t6_medium_silicate_slowsw_betaderivation.npz" , x = x[::10] , y = y[::10] , vx = vx[::10] , vy = vy[::10] , m = m[::10] , b = b[::10] , t = t[::10])

    # x , y , b , t = vals[: , 0] , vals[: , 1] , vals[: , 5] , vals[: , 6]
    # np.savez("Files/rk45_t6_large_carbon_slowsw_31.6AU.npz" , x = x[::10] , y = y[::10] , b = b[::10] , t = t[::10])
    
    plt.plot(t[::10] , np.sqrt(x**2+y**2)[::10])
    plt.show()
    
    