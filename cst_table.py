import numpy as np
from scipy.constants import *
from constants import *

"""for calculating epsilon for mass loss by selecting physically meaningful
mass values....."""

"""class describing particle's attributes for dust particle in bound orbit
around the Sun"""
class particle_scaled:
    """attributes: rho (float), assumed particle density
                   radius (float), set particle radii
                   
       properties: mass (float), particle mass in kg calculated from radius
                   and density
                   
                   area (float), cross section calcs based on radius
                   beta (float), beta calcs based on mass and initial cond
                   period (float), orbital period calcs based on beta"""
    
    def __init__(self , radius , beta):
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
    
    """beta calcs"""
    def beta(self):
        m = self.mass()
        m0 = m.flat[0]
        
        b = beta_const * self.b0 * (m0 / m)**(1 / 3) #scaled beta formula
        
        return b
    
    """period calcs"""
    def period(self):
        b = self.beta()
        
        T = period_cst * np.sqrt(1 / (1 - b)) #orbital period
        
        return T
        

if __name__ == "__main__":
    rad0 = np.linspace(500e-9 , 10e-6 , 10) #large initial particle radii, m
    
    particles = [particle_scaled(r , beta_0) for r in rad0] #particle attributes
    masses = np.array([p.mass() for p in particles]) #mass array
    m0 = masses[0]
    
    beta = np.array([0.45931934 , 0.14763836 , 0.08795477 , 0.06263446 , 0.04863381 ,
                     0.03974879 , 0.03360873 , 0.02911179 , 0.02567624 , 0.02296597])
    
    T = period_cst * np.sqrt(1 / (1 - beta))
    
    tot = mass_cst * T * m0**(-1 / 3)
    print(T / yr , tot)
    
    
    
    
    
    



