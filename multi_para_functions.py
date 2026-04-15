import numpy as np
from config import t5 , t6 , t7
from dust_properties import dust_properties

class multi_parameter():
    def __init__(self , t , particle):
        self.time = t
        self.B = particle.B
        self.epsilon = particle.epsilon
        self.delta = particle.delta

    def r(self):
        r0 = 1.0
        r10 = self.B / (3 * (1 - self.B)) * (self.time - np.sin(self.time))
        r01 = -self.B / (3 * (1 - self.B)) * (np.sin(self.time) + 5 * self.time)

        rtot = r0 + self.epsilon * r10 + self.delta * r01

        return rtot
    
    def theta(self):
        theta0 = self.time
        theta10 = self.B / (3 * (1 - self.B)) * (2 - 2 * np.cos(self.time) - self.time**2)
        theta01 = self.B / (3 * (1 - self.B)) * (7 * self.time**2 / 2 - 2 * np.cos(self.time) + 2)

        thetatot = theta0 + self.epsilon * theta10 + self.delta * theta01

        return thetatot
    
    def vr(self):
        vr0 = 0.0
        vr10 = self.B / (3 * (1 - self.B)) * (1 - np.cos(self.time))
        vr01 = -self.B / (3 * (1 - self.B)) * (np.cos(self.time) + 5)

        vrtot = vr0 + self.epsilon * vr10 + self.delta * vr01

        return vrtot

    def omega(self):
        omega0 = 1
        omega10 = 2 * self.B / (3 * (1 - self.B)) * (np.sin(self.time) - self.time)
        omega01 = self.B / (3 * (1 - self.B)) * (7 * self.time + 2 * np.sin(self.time))
        
        omegatot = omega0 + self.epsilon * omega10 + self.delta * omega01

        return omegatot
    
    def beta(self):
        beta0 = 1.0
        beta10 = self.time / 3
        beta01 = 0.0

        betatot = beta0 + self.epsilon * beta10 + self.delta * beta01

        return betatot

if __name__ == "__main__":
    par = dust_properties("carbon" , "slow" , "all" , "large")
    num = np.load("Files/rk45_t6_large_carbon_slowsw.npz")
    _ , _ , _ , _ , _ , b , t = [num[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]

    res = multi_parameter(t , par)

    rtot , thetatot , vrtot , omegatot , betatot = res.r() , res.theta() , res.vr() , res.omega() , res.beta()

    np.savez("Files/multiparam_t6_large_carbon_slowsw.npz" , omega = omegatot , r = rtot , theta = thetatot , vr = vrtot , c0 = 0.0 , b = betatot , t = t)