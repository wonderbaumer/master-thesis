from scipy.integrate import solve_ivp
import numpy as np
from tqdm import tqdm
from config import t6 , t7 , t8 , t9 , t10
from forces_scaled import tot_acc, sputtering, betahat 

#Source for pbar is Lima, 2020

"""ode function that is put into particle motion ivp solver, on form
dxdt, dydt, dvxdt, dvydt, dmdt"""
def pos_vel(t , init , pbar , state , particle_obj , massloss = True , drag = True):
    """input: t (float), time in s
              init (array), initial conditions x, y, vx, vy, m
              massloss, default: True, specify False if massloss is not considered
              
        returns: variable_list (list), 
        list of calculated variables for ivp solver, dxdt, dydt, dvxdt, dvydt, dmdt"""

    last_t , dt = state
    
    # let's subdivide t_span into 1000 parts
    # call update(n) here where n = (t - last_t) / dt

    n = int((t - last_t) / dt)
    pbar.update(n)
    ###
    # this we need to take into account that n is a rounded number.
    state[0] = last_t + dt * n
    
    x , y , vx , vy , m  = init #initial conditions unpacked

    dmdt = sputtering(m , particle_obj.epsilon , x , y) if massloss else 0.0 #specify change in mass if massloss is considered

    ax , ay = tot_acc(x , y , vx , vy , m , particle_obj , drag) #acceleration calcs

    variable_list = np.array([vx , vy , ax , ay , dmdt]) #variables for ivp solver
    
    return variable_list

"""Stopping condition if size outside interpolation range"""
def r_out_of_range_event(t , init , pbar , state , particle_obj , massloss , drag):
    _ , _ , _ , _ , m = init

    lower_lim = 10**(-23)
    m_physical = m * particle_obj.m0
    
    return m_physical - lower_lim

r_out_of_range_event.terminal = True
r_out_of_range_event.direction = -1

"""Stopping condition if distance <0.1"""
def orbital_radius_event(t , init , pbar , state , particle_obj , massloss , drag):
    x , y , _ , _ , _ = init

    r_orbit = np.sqrt(x**2 + y**2)

    return r_orbit - 1e-2

orbital_radius_event.terminal = True
orbital_radius_event.direction = -1
    
"""particle motion solved using non-stiff solver in scipys solve_ivp"""
def particle_motion(fun , t_span , y0 , method , state , particle_obj , massloss = True , drag = True):
    """input: fun (func), system of ODEs to solve
              t_span (tuple), time start and end of integration
              y0 (array), initial conditions for pos and vel
              method (string), solver to use, options RK45, RK23, DOP853
              state (array), updating time for pbar
              particle_obj (instance), containing particle information
              massloss (bool), True or False depending on whether massloss is considered
              
       returns: sol (array_like), x, y, vx and vy, m at time t"""
    
    with tqdm(total = 1000) as pbar:
        sol = solve_ivp(fun , t_span , y0 , method = method ,
                      args = (pbar , state , particle_obj , massloss , drag) , rtol = 1e-6 
                      , atol = 1e-9 , events = [r_out_of_range_event , orbital_radius_event]) #solving diff eq using solve_ivp, tight tolerances

    return sol

"""function that creates array of variables from solution object"""
def arr_variables(sol , particle_obj , massloss = True):
    """input: sol (array_like), solution object from solve_ivp
              particle_obj (instance), containing particle information
              
       
       returns: new_arr (array), array containing x , y , vx , vy , m , beta values"""
    
    x , y , vx , vy , m = sol.y #unpacking solution object
    b = betahat(m , particle_obj)

    new_arr = np.column_stack((x , y , vx , vy , m , b , sol.t)) #creating new array with all variables

    return new_arr

    