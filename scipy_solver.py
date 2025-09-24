from scipy.integrate import solve_ivp
import numpy as np
from polar_to_cart import polar_to_cartesian
from forces import *

"""ode function that is put into particle motion ivp solver, on form
dxdt, dydt"""
def pos_vel(t , init , r_par):
    """input: t (float), time in s
              init (array), initial conditions position and vel in x and y
              r_par (float), radius of particle
              acc (func), acceleration function of particle
              
        returns: variable_list (list), 
        list of calculated vel and acc in x and y"""
        
    x , y , vx , vy = init #initial conditions
    ax , ay = tot_acc(x , y , r_par) #acceleration unpacking
    
    variable_list = [vx , vy , ax , ay] #calculated vel and acc
    
    return variable_list
    
"""particle motion solved using non-stiff solver in scipys solve_ivp"""
def particle_motion(fun , t_span , y0):
    """input: fun (func), system of ODEs to solve
              t_span (tuple), time start and end of integration
              y0 (array), initial conditions for pos and vel
              
       returns: sol (array_like), x ,y ,vx and vy at time t"""
    
    sol = solve_ivp(fun , t_span , y0 , method = "DOP853" , args=(r_par,) ,
                    rtol=1e-9 , atol=1e-12)
    
    return sol

if __name__ == "__main__":
    theta0 = 0 #initial angle in rad, initial position along horizontal
    v0r = 0 #initial radial vel in m/s
    v0theta = 29.78e3 #initial angular vel in m/s
    
    init_polar = np.array([r0 , theta0 , v0r , v0theta]) #initial values array
    init_cartesian = polar_to_cartesian(init_polar) #initial values to cartesian
    
    t0 = 0 #initial time in s
    t_tot = 3.16e10 #total time in s
    t_span = (t0 , t_tot) #tuple of start and end time

    pos_and_vel = particle_motion(pos_vel , t_span , init_cartesian[0])
    
    
   