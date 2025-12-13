import numpy as np
import sys
sys.path.insert(1, 'C:/Users/Cecilie.Bamer/Documents/Project-paper/')
from forces import *
from polar_to_cart import *
from scipy.constants import *
from config import *
from leapfrog import *
from scipy_solver import *

sim_label = {t1 : "t1" , t2 : "t2" , t3 : "t3" , t4 : "t4"} #for filenames

"""class calculating and plotting the particle's orbit and energy, 
includes forces, mass loss processes, user-specified solver"""
class particle():
    """attributes:
       
       init_cond (array), initial radial position, theta angular position and
                          radial velocity 
                          
       sim_time (tuple), consisting of dt and t_tot, where dt is timestep used
                         in solver, t_tot is total time used in solver
                         dt should be smaller than or equal to T/20
        
       solver (string), user-specified solver for trajectory. 
                        Default:LEAPFROG, else RK45, RK23 or DOP853
                   
       m (float), optional: default: m_par, user can specify different mass than "standard" used in sims

       massloss (boolean), optional. default:true as massloss considered, false if not
                   
       methods: 
       
       initial_vel(), use initial radial distance and specified mass to calculate initial angular velocity,
       adds it to the initial values array, returns the total initial values
       
       init_conds_cart(), converting initial values to cartesian coordinates, returns cartesian coordinates
       
       beta(x , y , m), calculates ratio between pressure radiation and gravitational force, based on position
                        and mass, returns the beta value
       
       acceleration(x , y , m), calculates total acceleration of particle based on position and mass,
                                returns acceleration value
        
       ode_vars(t , init), provides parameters for scipy ivp solver in form vx, vy, ax, ay, dm
       
       pos_vel_calcs(), numerical solution for particle parameters based on input of solver type and 
                        if massloss or not"""
    
    def __init__(self , init_cond , sim_time , solver = "LEAPFROG" , m = m_par , massloss = True):
        
        self.m = m
        self.init_cond = init_cond
        self.solver = solver
        self.sim_time = sim_time
        self.massloss = massloss
        self.sim_label = sim_label[sim_time]
    
    """calculates necessary initial angular velocity for a particle at a 
    specified distance from larger object, for the particle to stay in orbit"""
    def initial_vel(self):
        
        r0 , theta0 , vr0 = self.init_cond
        
        vtheta0 = np.sqrt(G * m_s * (1 - beta0) / r0)
        
        """calculates and returns full set of polar initial conditions in 
        sequence radial distance, angular distance, radial velocity, angular
        velocity"""
        self.tot_init = np.append(self.init_cond , vtheta0)
        
        return self.tot_init
    
    """converts the array of initial conditions from polar to cartesian"""
    def init_conds_cart(self):
        
        init_conds_polar = self.initial_vel()
        
        init_conds_cart = polar_to_cartesian(init_conds_polar)
        
        return init_conds_cart 
    
    """beta calcs based on mass of particle and forces affecting it"""
    def beta(self , x , y , m):
        
        gx , gy = gravity(x , y)
        g_abs = np.sqrt(gx**2 + gy**2)
        
        Fr_abs = pressure_radial(x , y , m)
        Fr_acc = Fr_abs / m
        
        beta_val = Fr_acc / g_abs
        
        return beta_val
    
    """function calculating acceleration based on forces.
    if user specifies beta=0, acceleration is only based on Fg, if beta not
    specified then acceleration is calculated from Fg and Fr"""
    def acceleration(self , x , y , m):
        """input: x (float), x position
                  y (float), y position
           
            returns: a (tuple), acceleration in x and y direction """
            
        a = tot_acc(x , y , m) #Fg and Fr
        
        return a
    
    """function that calculates velocity and acceleration as function of
    time, for scipy solvers"""
    def ode_vars(self , t , init):
        """input: t (float), time in seconds
                  init (list), initial values for position and velocity
                  m (float), mass in kg
                  
          returns: ode_variables (list), dx, dy, ddx, ddy"""
        
        x , y_pos , vx , vy , m = init #initial value unpacking
        dm = sputtering(m) if self.massloss else 0.0 #define mass change per time parameter, 0 if no massloss
        ax , ay = self.acceleration(x , y_pos , m) #acceleration, ax, ay unpacking
        
        ode_variables = [vx , vy , ax , ay , dm] #values on form vx, vy, ax, ay, dm
        
        return ode_variables
    
    """calculates position, velocity and other parameters using different solvers"""
    def pos_vel_calcs(self):
        
        m = self.m #mass
        
        initial_vals = self.init_conds_cart() #initial values for leapfrog
        y0 = np.append(initial_vals , m) #initial values for scipy ivp solver
        
        dt , t_tot = self.sim_time #dt and t_tot unpacking
        
        t_span = (0 , t_tot) #time for simulations
        t_eval = np.arange(0 , t_tot , dt) #setting number of timesteps scipy solver
        
        if self.solver == "LEAPFROG" and self.massloss == True:
            pos_and_vel = leapfrog_algorithm(initial_vals , self.acceleration
                     , dt , t_tot , sputtering) #leapfroging using initial cond
            
        elif self.solver == "LEAPFROG" and self.massloss == False:
            pos_and_vel = leapfrog_algorithm(initial_vals , self.acceleration
                     , dt , t_tot ) #leapfroging using initial cond
        else:
            pos_and_vel = particle_motion(self.ode_vars , t_span , 
                                          y0 , self.solver , t_eval) #specified scipy solver
            
            x , y , vx , vy , m = pos_and_vel.y[0 , :] , pos_and_vel.y[1 , :] , pos_and_vel.y[2 , :] , pos_and_vel.y[3 , :] , pos_and_vel.y[4 , :]
            b_vals = np.array([self.beta(x[i] , y[i] , m[i]) for i in range(len(x))]) #beta calcs
            pos_and_vel = [x , y , vx , vy , m , b_vals] #new array of values

        return pos_and_vel
    
    def save_to_file(self):
        vals = self.pos_vel_calcs()

        if self.solver == "LEAPFROG":
            pos , b_vals = vals
            x , y , vx , vy , m = pos[: , 0] , pos[: , 1] , pos[: , 2] , pos[: , 3] , pos[: , 4]
            np.savez(f"C:/Users/Cecilie.Bamer/Documents/Project-paper/Files/{self.solver}_{self.sim_label}_massloss{self.massloss}.npz" , x = x , y = y , vx = vx , vy = vy , m = m , b = b_vals)

        else:
           x , y , vx , vy , m , b_vals = vals
           np.savez(f"C:/Users/Cecilie.Bamer/Documents/Project-paper/Files/{self.solver}_{self.sim_label}_massloss{self.massloss}.npz" , x = x , y = y , vx = vx , vy = vy , m = m , b = b_vals)

if __name__ == "__main__":

    p = particle(init_polar , t4 , "RK45")
    #p_rk = p.save_to_file()
    
    
    
    
    
        