from scipy.integrate import solve_ivp
import numpy as np
from tqdm import tqdm
from AstronomicalSilicate_modified import sputter
from config import t5 , t6 , t7 , init_cart_scaled , mhat0 , eps
from forces_scaled import tot_acc, sputtering, betahat 

#Source for pbar is Lima, 2020

"""ode function that is put into particle motion ivp solver, on form
dxdt, dydt, dvxdt, dvydt, dmdt"""
def pos_vel(t , init , pbar , state , epsilon , massloss = True):
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
    
    x , y , vx , vy , m = init #initial conditions unpacked
    dmdt = sputtering(m , epsilon) if massloss else 0.0 #specify change in mass if massloss is considered
    ax , ay = tot_acc(x , y , m) #acceleration calcs

    variable_list = np.array([vx , vy , ax , ay , dmdt]) #variables for ivp solver
    
    return variable_list
    
"""particle motion solved using non-stiff solver in scipys solve_ivp"""
def particle_motion(fun , t_span , y0 , method , t_eval , state , epsilon , massloss = True):
    """input: fun (func), system of ODEs to solve
              t_span (tuple), time start and end of integration
              y0 (array), initial conditions for pos and vel
              method (string), solver to use, options RK45, RK23, DOP853
              
       returns: sol (array_like) , x , y , vx and vy , m at time t"""
    
    with tqdm(total = 1000) as pbar:
        sol = solve_ivp(fun , t_span , y0 , method = method ,
                    t_eval = t_eval , args = (pbar , state , epsilon , massloss , ) , rtol=1e-9 , atol=1e-12) #solving diff eq using solve_ivp, tight tolerances
    
    return sol

"""function that creates array of variables from solution object, same structure as leapfrog output"""
def arr_variables(sol):
    """input: sol (array_like), solution object from solve_ivp
       
       returns: new_arr (array), array containing x , y , vx , vy , m , beta values"""
    
    x , y , vx , vy , m = sol.y #unpacking solution object
    beta = betahat(m) #calculating beta values from mass array

    new_arr = np.column_stack((x , y , vx , vy , m , beta)) #creating new array with all variables

    return new_arr


if __name__ == "__main__":
    dt , t_tot = t6
    t = np.arange(0 , t_tot , dt)
    state = [0 , t_tot / 1000]
    y0 = np.append(init_cart_scaled , mhat0)
    epsilon = eps("slow" , "all")


    sol = particle_motion(pos_vel , (0 , t_tot) , y0 , "RK45" , t , state , epsilon , massloss = True)
    lf_vals  = arr_variables(sol)
    x , y , vx , vy , m , b = lf_vals[: , 0] , lf_vals[: , 1] , lf_vals[: , 2] , lf_vals[: , 3] , lf_vals[: , 4] , lf_vals[: , 5]
   
    
    