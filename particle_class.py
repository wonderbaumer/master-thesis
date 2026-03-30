import numpy as np
from config import mhat0 , init_cart_scaled , init_cart , t5 , t6 , t7 , eps , r1 , r2 , r3 , rho
from leapfrog import leapfrog_algorithm
from scipy_solver import particle_motion , pos_vel , arr_variables
from forces_scaled import tot_acc, sputtering
from constants import m_s, au, c
from scipy.constants import G
from polar_to_cart import polar_to_cartesian


init_vals = {"large":{
            "r": 1.54079 * 10**(-6),
            "B": {"silicate": 0.1235 , "carbon": 0.2646}
            },

            "medium":{
            "r": 0.17508 * 10**(-6),
            "B": {"silicate": 0.8560 , "carbon": 3.0589}
            },

            "small":{
            "r": 0.04259 * 10**(-6),
            "B":{"silicate": 0.2098 , "carbon": 1.6179}     
            }

        }

"""class solving scaled equations of motion for a particle using user-specified numerical solver"""
class particle():
    """attributes:                
       sim_time (tuple), consisting of dt and t_tot
        
       solver (string), user-specified solver for equation of motion. 
                        Default:LEAPFROG, else RK45, RK23 or DOP853

       massloss (boolean), optional. default:True if massloss considered, False if not
                   
       methods: 
       
       pos_vel_calcs(), numerical solution for particle parameters based on input of solver type and 
                        if massloss=True or =False"""
    
    def __init__(self , sim_time , size ,  solver = "LEAPFROG" , massloss = True , material = "Silicate" , sw = "slow"):
        self.solver = solver
        self.sim_time = sim_time
        self.massloss = massloss
        self.material = material
        self.size = size
        self.sw = sw
    
        self.r , self.B = init_vals[self.size][self.material]
        self.epsilon = eps(self.m() , material , sw , "all")
    
    def m(self):
        m = 4 / 3 * np.pi * rho * (self.r)**3

        return m
    
    def init_params(self):
        R = 1 * au #initial radial position
        V = np.sqrt((G * m_s * (1 - self.B)) / R) #initial angular velocity, scaled formula
        self.T = round(np.sqrt(R**3 / (G * m_s * (1 - self.B)))) #initial period, scaled formula
        self.K = 1

        """scaled initial parameters"""
        rhat0 = R / R #initial scaled radial position
        thetahat0 = 0 #initial scaled angular position
        vrhat0 = 0  #initial scaled radial velocity
        omegahat0 = V / V #initial scaled angular velocity
        self.betahat0 = self.B / self.B #initial scaled beta
        self.mhat0 = self.m / self.m #initial scaled mass
        self.delta = V / c

        init_polar = np.array([R , 0 , 0 , V]) #initial polar coords
        self.init_cart = polar_to_cartesian(init_polar) #initial cart coords

        init_polar_scaled = np.array([rhat0 , thetahat0 , vrhat0 , omegahat0]) #initial scaled polar coords
        init_cart_scaled = polar_to_cartesian(init_polar_scaled) #initial scaled cart coords
        
        return init_cart_scaled

    """calculates position, velocity and other parameters using different solvers"""
    def pos_vel_calcs(self):
        y0 = np.append(self.init_cart_scaled, self.mhat0) #initial values for scipy ivp solver
        dt , t_tot = self.sim_time #dt and t_tot unpacking
        
        t_span = (0 , t_tot) #time for simulations
        t_eval = np.arange(0 , t_tot , dt) #setting number of timesteps scipy solver

        state = [0 , t_tot / 1000]
        
        if self.solver == "LEAPFROG" and self.massloss == True:
            pos_and_vel1 = leapfrog_algorithm(init_cart_scaled , tot_acc
                     , self.sim_time , self.epsilon , sputtering , self.material) #leapfroging using initial cond
            
        elif self.solver == "LEAPFROG" and self.massloss == False:
            pos_and_vel1 = leapfrog_algorithm(init_cart_scaled , tot_acc
                     , self.sim_time , material = self.material) #leapfroging using initial cond
            
        elif self.solver in ["RK45" , "RK23" , "DOP853"] and self.massloss == True: 
            pos_and_vel = particle_motion(pos_vel , t_span , 
                                          y0 , self.solver , t_eval , state , self.epsilon ,  massloss = True) #specified scipy solver
            
            pos_and_vel1 = arr_variables(pos_and_vel , self.material) #array of variables
        
        elif self.solver in ["RK45" , "RK23" , "DOP853"] and self.massloss == False: 
            pos_and_vel = particle_motion(pos_vel , t_span , y0 , self.solver , t_eval , state , self.epsilon, massloss = False) #specified scipy solver
            
            pos_and_vel1 = arr_variables(pos_and_vel , self.material) #array of variables

        return pos_and_vel1

if __name__ == "__main__":
    
    epsilon = eps(material = "carbon")
    p = particle(t5 , epsilon , "RK45" , massloss = True , material = "Carbon")
    vals = p.pos_vel_calcs()
    x , y , vx , vy , m , b , t = vals[: , 0] , vals[: , 1] , vals[: , 2] , vals[: , 3] , vals[: , 4] , vals[: , 5] , vals[: , 6]

    np.savez("Files/rk45_t5_carbon_slowsw_realbeta_test.npz" , x = x , y = y , vx = vx , vy = vy , m = m , b = b , t = t)
    