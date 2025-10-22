import numpy as np
from scipy.constants import *
from constants import *
import matplotlib.pyplot as plt

"""class describing particle's attributes for dust particle in bound orbit
around the Sun"""
class particle_scaled:
    """attributes: 
       rho (float), assumed particle density from constants
       radius (array), set particle radii
       beta_0 (array), beta_0 values
                   
       properties: 
       mass (): returns array of particle mass
       area (): returns array cross section 
       epsilon_mass (): returns array of epsilon values"""
    
    def __init__(self , radius , beta_0):
        self.r = radius
        self.rho = rho
        self.b0 = beta_0
    
    """mass calcs based on constant density and perfect sphere"""
    def mass(self):
        m = 4 * np.pi * self.rho * self.r**3 / 3 #mass in kg
        
        return m 
    
    """area calcs of sphere, cross section"""
    def area(self):
        A = np.pi * self.r**2
        
        return A
    
    """epsilon calcs"""
    def epsilon_mass(self):
        area_fraction = 3 / (4 * np.pi * self.rho)
        
        m = self.mass()
        self.m0 = m[0]
        
        T = period_cst * np.sqrt(1 / (1 - self.b0)) #orbital period
        
        eps = fsw * Ytot * mA * np.pi * T * area_fraction**(2 / 3) * self.m0**(-1 / 3)
        
        return eps
    
    """plotting function"""
    def plot(self):
        x = self.b0
        y = self.epsilon_mass()
        
        plt.plot(x , y)
        
        plt.xlabel("Beta_0")
        plt.ylabel("Epsilon")
        plt.title("Eps as function of beta_0")
        
        plt.show()
        
if __name__ == "__main__":
    rad0 = np.linspace(500e-9 , 10e-6 , 10) #large initial particle radii, m
    beta_0 = np.array([0.45931934 , 0.14763836 , 0.08795477 , 0.06263446 , 0.04863381 ,
                     0.03974879 , 0.03360873 , 0.02911179 , 0.02567624 , 0.02296597])
    
    p = particle_scaled(rad0 , beta_0)
    mass = p.mass()
    epsilon = p.epsilon_mass()
    print(beta_0)
    
    
    
    
    
    
    
    
    
    



