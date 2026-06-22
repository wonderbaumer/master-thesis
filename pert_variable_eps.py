from dust_properties import dust_properties
from scipy.integrate import solve_ivp
from config import t6 , t7 , t8 , t9 , t10 , yr , lower_m_sil , lower_m_car
import numpy as np
import matplotlib.pyplot as plt
import math


"""This code solves the zeroth order perturbed expressions (multiscale) numerically"""


class PerturbationSolver:
    """
    Solves zeroth-order multiscale perturbation equations.
    """

    def __init__(self, sim_time ,particle):

        self.time = sim_time
        self.par = particle


        self.B = particle.B
        self.K = particle.K

        self.y0 = [1.0, 1.0]

    # -----------------------------
    # ODE system
    # -----------------------------
    def rhs(self, time, y):
        beta0, C0 = y
        
        # m0 = 1.0
        # beta0 = m0**(-1 / 3)

        r0 = C0**2 * (1 - self.B) / (1 - self.B * beta0)
        
        dC0_dt = (-self.B * self.K * beta0 * (1 - beta0 * self.B)**2 / (1 - self.B)**3 * C0**(-3))

        dbeta0_dt = beta0**2 / 3 * r0**(-2)

        

        return [dbeta0_dt, dC0_dt]

    # -----------------------------
    # EVENT FACTORIES (IMPORTANT FIX)
    # -----------------------------
    def r_out_of_range_event(self):

        def event(time, y):
            beta0, C0 = y
            return None #(1 - self.B * beta0) - 1e-1

        event.terminal = True
        event.direction = -1
        return event

    def orbital_radius_event(self):

        def event(t, y):
            beta0, C0 = y
            
            r0 = C0**2 * (1 - self.B) / (1 - self.B * beta0)
            return r0 - 1e-1

        event.terminal = True
        event.direction = -1
        return event

    # -----------------------------
    # POST-PROCESSING
    # -----------------------------
    def arr_variables(self, sol):

        beta0, C0 = sol.y
        
        r0 = C0**2 * (1 - self.B) / (1 - self.B * beta0)
        
        omega0 = (
            C0**(-3)
            * ((1 - self.B) / (1 - self.B * beta0))**(-2)
        )

        # mpert = beta0**(-3)

        dbeta0_dt1 = beta0**2 / 3 * r0**(-2)


        return r0, omega0, beta0, sol.t , C0

    # -----------------------------
    # SOLVER
    # -----------------------------
    def solve(self):

        
        # print(teval , len(teval))
        # self.time = np.arange(0 , 30000 , 3.16e-3)
        teval = self.time
        sol = solve_ivp(
            self.rhs,
            (self.time[0], self.time[-1]),
            self.y0,
            # t_eval= self.time,
            method="RK45",
            rtol=1e-9,
            atol=1e-12,
            events=[
                #self.r_out_of_range_event(),
                self.orbital_radius_event(),
            ],
        )

        stopping_reason = None
        stopping_events = [
            #"particle outside interpolation range",
            "particle impacted Sun",
        ]

        for i, te in enumerate(sol.t_events):
            if len(te) > 0:
                stopping_reason = stopping_events[i]
                break

        print(stopping_reason)

        return self.arr_variables(sol)

if __name__== "__main__":
    material = "carbon"
    par = dust_properties(material , "fast" , init_dist = 1 , size = "A")
    res = np.load("Files/rk45_t8_C_silicate_slowsw.npz")
    x , y , _ , _ , m , b , t  = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    dt , t_tot = t10
    tarr = np.arange(0 , t_tot , dt)
    pert_solver = PerturbationSolver(
        tarr * par.epsilon,
        particle=par)

    r0, omega0, beta0, tpert , C0 = pert_solver.solve()

    # k = beta0 * (1 - par.B) / (6 * (1 - beta0 * par.B))
    # print(k[:-1])
    # plt.plot(tpert / par.epsilon , r0)
    
    # plt.show()
    print(r0[-1] , beta0[-1]**(-3) * par.m0 , tpert[-1] / par.epsilon * par.T / yr)
    
    
    