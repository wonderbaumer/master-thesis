from forces import *
from polar_to_cart import polar_to_cartesian
import numpy as np
from scipy.constants import *

"""class calculating and plotting the particle's orbit and energy, 
includes forces, mass loss processes, user-specified solver"""
class particle():
    """attributes:
       radius (float), user-specified radius of particle [m]
       
       init_cond (array), initial radial position, theta angular position and
                          radial velocity 
       
       beta (int), optional. Default:none and Fr and Fg calculated, 
                   else beta=0 and only Fg included in trajectory calcs
                   
       solver (string), user-specified solver for trajectory. 
                        Default:LEAPFROG, else RK45, RK23 or DOP853
        
       properties: """
    
    def __init__(self , radius , init_cond , solver = "LEAPFROG" , beta = None):
        self.r = radius
        self.cond = init_cond
        self.b = beta 
        self.sol = solver
    
    """calculates necessary initial angular velocity for a particle at a 
    specified distance from larger object, for the particle to stay in orbit"""
    
    def initial_vel(self):
        r0 , theta0 , vr0 = self.cond
        
        vtheta0 = np.sqrt(G * m_s / r0)
        
        self.tot_init = np.append(self.cond , vtheta0)
        
        return self.tot_init
    
    """mass calcs based on constant density and perfect sphere"""
    def mass(self):
        m = 4 * np.pi * rho * self.r**3 / 3 #mass in kg
        
        return m 
    
    """beta calcs based on mass of particle and forces affecting it"""
    def betas(self , x , y):
        
        m_par = self.mass()
        
        gx , gy = gravity(x , y)
        g_abs = np.sqrt(gx**2 + gy**2)
        
        Fr_abs = pressure_radial(x , y , self.r)
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
            
        beta_val = self.b if self.b is not None else self.betas(x , y)
            
        if beta_val != 0:
            a = tot_acc(x , y , beta_val) #Fg and Fr
        else:
            a = gravity(x , y) #Only Fg
        
        return a
    
    
            
            
            
    
    
    
    
        