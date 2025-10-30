from forces import *
from polar_to_cart import polar_to_cartesian
import numpy as np
from scipy.constants import *
from leapfrog import leapfrog_algorithm

"""class calculating and plotting the particle's orbit and energy, 
includes forces, mass loss processes, user-specified solver"""
class particle():
    """attributes:
       radius (float), user-specified radius of particle [m]
       
       init_cond (array), initial radial position, theta angular position and
                          radial velocity 
                          
       sim_time (tuple), consisting of dt and t_tot, where dt is timestep used
                         in solver, t_tot is total time used in solver
                         dt should be smaller than or equal to T/20
        
       solver (string), user-specified solver for trajectory. 
                        Default:LEAPFROG, else RK45, RK23 or DOP853
                        
       beta (int), optional. Default:none and Fr and Fg calculated, 
                   else beta=0 and only Fg included in trajectory calcs
                   
       
        
       properties: """
    
    def __init__(self , radius , init_cond , sim_time , solver = "LEAPFROG" , beta = None):
        self.radius = radius
        self.init_cond = init_cond
        self.beta_user = beta 
        self.solver = solver
        self.sim_time = sim_time
    
    """calculates necessary initial angular velocity for a particle at a 
    specified distance from larger object, for the particle to stay in orbit"""
    def initial_vel(self):
        r0 , theta0 , vr0 = self.init_cond
        
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
        
    
    """mass calcs based on constant density and perfect sphere"""
    def mass(self):
        m = 4 * np.pi * rho * self.radius**3 / 3 #mass in kg
        
        return m 
    
    """beta calcs based on mass of particle and forces affecting it"""
    def beta(self , x , y):
        
        m_par = self.mass()
        
        gx , gy = gravity(x , y)
        g_abs = np.sqrt(gx**2 + gy**2)
        
        Fr_abs = pressure_radial(x , y , self.radius)
        Fr_acc = Fr_abs / m_par
        
        beta_val = Fr_acc / g_abs
        
        return beta_val
    
    """function calculating acceleration based on forces.
    if user specifies beta=0, acceleration is only based on Fg, if beta not
    specified then acceleration is calculated from Fg and Fr"""
    def acceleration(self , x , y):
        """input: x (float), x position
                  y (float), y position
           
            returns: a (tuple), acceleration in x and y direction """
         
        beta_val = self.beta_user if self.beta_user is not None else self.beta(x , y)
            
        if beta_val == 0:
            a = gravity(x , y) #Only Fg
            
        else:
            a = tot_acc(x , y , beta_val) #Fg and Fr
        
        return a
    
    """calculates """
    def pos_vel_calcs(self):
        initial_vals = self.init_conds_cart()
        x , y , vx , vy = initial_vals
        
        acc_func = self.acceleration(x , y)
        dt , t_tot = self.sim_time
        
        
        
        if self.solver == "LEAPFROG":
            pos_and_vel = leapfrog_algorithm(init_cartesian , tot_acc , dt , t_tot) #leapfroging using initial cond
        
        else:
            pos_and_vel = particle_motion( , self.sim_time , 
                                          init_cartesian[0])
        
        return pos_and_vel
    
if __name__ == "__main__":
    init_polar = np.array([r0 , theta0 , v0r]) #initial values array
            
            
            
    
    
    
    
        