from scipy.integrate import solve_ivp
import numpy as np
from polar_to_cart import polar_to_cartesian
from forces import *
from constants import *
import matplotlib.pyplot as plt

"""ode function that is put into particle motion ivp solver, on form
dxdt, dydt"""
def pos_vel(t , init):
    """input: t (float), time in s
              init (array), initial conditions position and vel in x and y
              r_par (float), radius of particle
              acc (func), acceleration function of particle
              
        returns: variable_list (list), 
        list of calculated vel and acc in x and y"""
        
    x , y , vx , vy , m = init #initial conditions
    
    dm = sputtering(m) #mass change
    
    ax , ay = tot_acc(x , y , m) #acceleration unpacking
    print(ax)
    
    variable_list = np.array([vx , vy , ax , ay , dm]) #calculated vel and acc

    return variable_list
    
"""particle motion solved using non-stiff solver in scipys solve_ivp"""
def particle_motion(fun , t_span , y0 , method , t_eval):
    """input: fun (func), system of ODEs to solve
              t_span (tuple), time start and end of integration
              y0 (array), initial conditions for pos and vel
              method (string), solver to use, options RK45, RK23, DOP853
              
       returns: sol (array_like), x ,y ,vx and vy,m at time t"""
    
    
    sol = solve_ivp(fun , t_span , y0 , method = method ,
                    t_eval = t_eval , rtol=1e-9 , atol=1e-12)
    
    return sol

if __name__ == "__main__":
    v0theta = 2.19013101e+04 #initial angular vel in m/s

    init_polar = np.array([r0 , theta0 , v0r , v0theta]) #initial values array
    init_cartesian = polar_to_cartesian(init_polar) #initial values to cartesian
    init_cartesian = np.append(init_cartesian , m_par)

    t0 = 0 #initial time in s
    t_tot = 3.16e10 #total time in s
    t_span = (t0 , t_tot) #tuple of start and end time
    t_eval = np.linspace(t0 , t_tot , 100000)

    pos_and_vel = particle_motion(pos_vel , t_span , init_cartesian , "DOP853" , t_eval)

    plt.plot(pos_and_vel.y[0] , pos_and_vel.y[1])
    #plt.show()
    
    
   