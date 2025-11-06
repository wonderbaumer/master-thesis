from forces import *
from polar_to_cart import polar_to_cartesian
import numpy as np
from scipy.constants import *
from leapfrog import leapfrog_algorithm
import matplotlib.pyplot as plt
from scipy_solver import *

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
                        
       beta (int), optional. Default:none and Fr and Fg calculated, 
                   else beta=0 and only Fg included in trajectory calcs
                   
      mass (float), user-specified initial mass of particle [kg]
                   
       
        
       properties: """
    
    def __init__(self , init_cond , sim_time , solver = "LEAPFROG" , beta = None 
                 , m = m_par):
        
        self.m = m
        self.init_cond = init_cond
        self.beta_user = beta 
        self.solver = solver
        self.sim_time = sim_time
    
    """calculates necessary initial angular velocity for a particle at a 
    specified distance from larger object, for the particle to stay in orbit"""
    def initial_vel(self):
        
        r0 , theta0 , vr0 = self.init_cond
        
        if beta is not None:
            vtheta0 = np.sqrt(G * m_s * (1 - beta_0) / r0)
        
        else:
            vtheta0 = np.sqrt(G * m_s / r0)
        
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
         
        beta_val = self.beta_user if self.beta_user is not None else self.beta(x , y , m)
            
        if beta_val == 0:
            a = gravity(x , y) #Only Fg
            
        else:
            a = tot_acc(x , y , m) #Fg and Fr
        
        return a
    
    """function that calculates velocity and acceleration as function of
    time, for scipy solvers"""
    def ode_vars(self , t , init):
        """input: t (float), time
                  init (list), initial values for position and velocity
                  m (float), mass in kg
                  
          returns: ode_variables (list), dx,dy,ddx,ddy"""
        
        x , y_pos , vx , vy , m = init
        dm = sputtering(m)
        ax , ay = self.acceleration(x , y_pos , m)
        
        ode_variables = [vx , vy , ax , ay , dm]
        
        return ode_variables
    
    """calculates position, velocity and other parameters using different solvers"""
    def pos_vel_calcs(self):
        
        m = self.m #mass
        
        initial_vals = self.init_conds_cart() #initial values
        y0 = np.append(initial_vals , m)
        
        dt , t_tot = self.sim_time
        
        t_span = (0 , t_tot)
        t_eval = np.linspace(0 , 3.16e10 , 100000) #forcing timesteps scipy solver
        
        if self.solver == "LEAPFROG":
            pos_and_vel = leapfrog_algorithm(initial_vals , self.acceleration
                    , sputtering , dt , t_tot) #leapfroging using initial cond
        
        else:
            pos_and_vel = particle_motion(self.ode_vars , t_span , 
                                          y0 , self.solver , t_eval)
        
        return pos_and_vel

sim_time = (3.16e5 , 3.16e10)    

p = particle(init_polar , sim_time , "RK45")


       
    
    
    
    
        