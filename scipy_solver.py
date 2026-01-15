from scipy.integrate import solve_ivp
import numpy as np
from forces import *
from forces_scaled import tot_acc, sputtering, betahat


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

if __name__ == "__main__":
    dt , t_tot = t6
    t = np.arange(0 , t_tot , dt) / T
    m = m_par / m_par
    y0 = np.append(init_cart_scaled , m)

    sol = particle_motion(pos_vel , (0 , t_tot / T) , y0 , "RK45" , t)
    x = sol.y[0]   # array of x positions
    y = sol.y[1]  # array of y positions
    vx = sol.y[2]  # array of x velocities
    vy = sol.y[3]  # array of y velocities
    m = sol.y[4]   # array of mass
    b = beta(x , y , m)

    
    b = betahat(m)
    #np.savez(f"C:/Users/cecil/Documents/Project-paper/Files/rk45_t6_masslossFalse_scaledeqs.npz" , x = x , y = y , vx = vx , vy = vy , m = m , b = b)
    