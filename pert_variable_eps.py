from dust_properties import dust_properties
from scipy.integrate import solve_ivp
from config import t6 , t7 , t8 , t9 , t10 , yr , lower_m_sil , lower_m_car
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import time


"""This code solves the zeroth order perturbed expressions (multiscale) numerically"""

def rhs(t , y0 , pbar , state , particle_obj):

    last_t , dt = state
    
    # let's subdivide t_span into 1000 parts
    # call update(n) here where n = (t - last_t) / dt

    n = int((t - last_t) / dt)
    pbar.update(n)
    ###
    # this we need to take into account that n is a rounded number.
    state[0] = last_t + dt * n

    beta0, C0 = y0 #beta0 same as E(t1)
        
    r0 = C0**2 * (1 - particle_obj.B) / (1 - particle_obj.B * beta0)
        
    dC0_dt = (-particle_obj.B * particle_obj.K * beta0 * (1 - beta0 * particle_obj.B)**2 
              / (1 - particle_obj.B)**3 * C0**(-3))

    dbeta0_dt = beta0**2 / 3 * r0**(-2)

    return [dbeta0_dt, dC0_dt]

def orbital_radius_event(t , y0 , pbar , state , particle_obj):
    beta0, C0 = y0
            
    r0 = C0**2 * (1 - particle_obj.B) / (1 - particle_obj.B * beta0)
    return r0 - 1e-1

orbital_radius_event.terminal = True
orbital_radius_event.direction = -1       

def solve(func , t , y0 , state , particle_obj):

    with tqdm(total = 1000) as pbar:
        start = time.perf_counter()
        sol = solve_ivp(func , (t[0] , t[-1]) , y0 , "RK45" , 
                        args = (pbar , state , particle_obj) , rtol = 1e-9 ,  atol = 1e-12 , 
                        events = [orbital_radius_event])
        end = time.perf_counter()

    print(f"Solver runtime: {end - start:.6f} seconds")

    return sol

def arr_variables(sol , particle_obj):

    beta0 , C0 = sol.y
        
    r0 = C0**2 * (1 - particle_obj.B) / (1 - particle_obj.B * beta0)
        
    omega0 = (C0**(-3) * ((1 - particle_obj.B) / (1 - particle_obj.B * beta0))**(-2))

    return r0 , omega0 , beta0 , sol.t , C0

class runner_class():
    def __init__(self , par , sim_time):
        self.B = par.B
        self.K = par.K
        self.sim_time = sim_time
    
    def solver(self):
        y0 = [1.0 , 1.0]

        dt , t_tot = self.sim_time #dt and t_tot unpacking
        dt = dt  
        t_tot = t_tot 
        t_span = (0 , t_tot) #time for simulations
        state = [0 , t_tot / 1000] #for progress bar


        stopping_reason = None
        stopping_events = ["particle impacted Sun"]

        orbit = solve(rhs , t_span , y0 , state , self)

        for i , te in enumerate(orbit.t_events):
            if len(te) > 0:
                stopping_reason = stopping_events[i]
                break

        print(stopping_reason)

        tot_variables = arr_variables(orbit , self)

        return tot_variables

if __name__== "__main__":
    par = dust_properties("silicate" , "CME" , init_dist = 1 , size = "C")
    p = runner_class(par , t6)
    vals = p.solver()
    
    r0 , omega0 , beta0 , t , C0 = vals[0] , vals[1] , vals[2] , vals[3] , vals[4]
    Kcst = beta0 * (1 - par.B)**2 / (6 * (1 - beta0 * par.B))
    # plt.plot(t / par.epsilon , r0)
    # plt.show()
    print(par.K)
    