from scipy.integrate import solve_ivp
import numpy as np
from tqdm import tqdm
from config import t5 , t6 , t7  
from forces_scaled import tot_acc, sputtering, betahat 

#Source for pbar is Lima, 2020

"""ode function that is put into particle motion ivp solver, on form
dxdt, dydt, dvxdt, dvydt, dmdt"""
def pos_vel(t , init , pbar , state , epsilon , particle_obj , massloss = True):
    """input: t (float), time in s
              init (array), initial conditions x, y, vx, vy, m
              massloss, default: True, specify False if massloss is not considered
              
        returns: variable_list (list), 
        list of calculated variables for ivp solver, dxdt, dydt, dvxdt, dvydt, dmdt"""

    last_t, dt = state
    
    # let's subdivide t_span into 1000 parts
    # call update(n) here where n = (t - last_t) / dt

    n = int((t - last_t)/dt)
    pbar.update(n)

    # this we need to take into account that n is a rounded number.
    state[0] = last_t + dt * n
    
    x , y , vx , vy , m  = init #initial conditions unpacked

    dmdt = sputtering(m , epsilon , x , y) if massloss else 0.0 #specify change in mass if massloss is considered

    ax , ay = tot_acc(x , y , vx , vy , m , particle_obj) #acceleration calcs

    variable_list = np.array([vx , vy , ax , ay , dmdt]) #variables for ivp solver
    
    return variable_list

def r_out_of_range_event(t , init , pbar , state , epsilon , particle_obj , massloss):
    x , y , vx , vy , m = init

    size = m**(1/3)
    r_physical = size * particle_obj.r 
    
    return r_physical - 1e-9

r_out_of_range_event.terminal = True
r_out_of_range_event.direction = -1

def orbital_radius_event(t, init, pbar, state, epsilon, particle_obj, massloss):
    x , y , vx , vy , m = init

    r_orbit = np.sqrt(x**2 + y**2)

    return r_orbit - 1e-3


orbital_radius_event.terminal = True
orbital_radius_event.direction = -1
    
"""particle motion solved using non-stiff solver in scipys solve_ivp"""
def particle_motion(fun , t_span , y0 , method , t_eval , state , epsilon , particle_obj , massloss = True):
    """input: fun (func), system of ODEs to solve
              t_span (tuple), time start and end of integration
              y0 (array), initial conditions for pos and vel
              method (string), solver to use, options RK45, RK23, DOP853
              
       returns: sol (array_like) , x , y , vx and vy , m at time t"""
    
    with tqdm(total = 1000) as pbar:
        sol = solve_ivp(fun , t_span , y0 , method = method ,
                      args = (pbar , state , epsilon , particle_obj , massloss , ) , rtol=1e-9 , atol=1e-12 , events = [r_out_of_range_event , orbital_radius_event]) #solving diff eq using solve_ivp, tight tolerances

    
    return sol

"""function that creates array of variables from solution object, same structure as leapfrog output"""
def arr_variables(sol , particle_obj):
    """input: sol (array_like), solution object from solve_ivp
       
       returns: new_arr (array), array containing x , y , vx , vy , m , beta values"""
    
    x , y , vx , vy , m = sol.y #unpacking solution object
    b = betahat(m , particle_obj)

    new_arr = np.column_stack((x , y , vx , vy , m , b , sol.t)) #creating new array with all variables

    return new_arr


if __name__ == "__main__":
    dt , t_tot = t5
    # t = np.arange(0 , t_tot , dt)
    # state = [0 , t_tot / 1000]
    #y0 = np.append(init_cart_scaled , mhat0)
    # epsilon = eps()


    # sol = particle_motion(pos_vel , (0 , t_tot) , y0 , "RK45" , t , state , epsilon , massloss = False)
    # lf_vals  = arr_variables(sol)
    # x , y , vx , vy , m , b = lf_vals[: , 0] , lf_vals[: , 1] , lf_vals[: , 2] , lf_vals[: , 3] , lf_vals[: , 4] , lf_vals[: , 5]
    # print(epsilon)
    
    