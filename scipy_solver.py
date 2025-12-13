from scipy.integrate import solve_ivp
import numpy as np
import sys
sys.path.insert(1, 'C:/Users/Cecilie.Bamer/Documents/Project-paper/')
from forces import *


"""ode function that is put into particle motion ivp solver, on form
dxdt, dydt, dvxdt, dvydt, dmdt"""
def pos_vel(t , init , massloss = True):
    """input: t (float), time in s
              init (array), initial conditions x, y, vx, vy, m
              massloss, default: True, specify False if massloss is not considered
              
        returns: variable_list (list), 
        list of calculated variables for ivp solver, dxdt, dydt, dvxdt, dvydt, dmdt"""

    x , y , vx , vy , m = init #initial position, velocity (cartesian coordinates) and mass

    dmdt = sputtering(m) if massloss else 0.0 #specify change in mass if massloss is considered
    ax , ay = tot_acc(x , y , m) #acceleration calcs

    variable_list = np.array([vx , vy , ax , ay , dmdt]) #variables for ivp solver
    
    return variable_list
    
"""particle motion solved using non-stiff solver in scipys solve_ivp"""
def particle_motion(fun , t_span , y0 , method , t_eval):
    """input: fun (func), system of ODEs to solve
              t_span (tuple), time start and end of integration
              y0 (array), initial conditions for pos and vel
              method (string), solver to use, options RK45, RK23, DOP853
              
       returns: sol (array_like) , x , y , vx and vy , m at time t"""
    
    
    sol = solve_ivp(fun , t_span , y0 , method = method ,
                    t_eval = t_eval , rtol=1e-9 , atol=1e-12) #solving diff eq using solve_ivp, tight tolerances
    
    return sol


    
    
    
   