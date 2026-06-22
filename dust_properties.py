from config import (au , init_vals , rho_s , rho_c , yr , m_s , sil_size ,
                     sil_betaval , car_betaval_bound , car_size_bound , mA_S , mA_C 
                     , size_to_mass , car_betaval , car_size)
from sputtering_dict import sputter
from scipy.constants import G , c
import numpy as np

class dust_properties():
    """Calculating properties of particle given speicified material, 
       solar wind conditions, initial distance and initial size
    
    Attributes: material (string), "silicate" or "carbon", particle material
                sw (string), "fast", "slow" or "CME", solar wind conditions
                size (string), default: None, "A", "C", "F", "G", "H", "D", "E", initial particle size
                size_range (tuple), default: None, else tuple consisting of a set of size values, beta values
                init_dist (float), initial distance of particle in AU

                r (float), initial particle size, in m
                B (float), initial beta value, unitless
                m0 (float), initial particle mass, in kg
                delta (float), PR drag term, unitless
                K (float), PR drag to mass loss ratio, unitless

    Methods: calc_V(), calculates initial orbital velocity, ms^-1
             calc_T(), calculates initial orbital period, s
             sw_flux(), calculates solar wind flux, kgm^-3
             sputtering_yield(), calculates total sputtering yield, atoms/incident ions
             sputtering_lifetime(), calculates sputtering lifetime
             eps(), calculates epsilon, mass loss rate, unitless
                """
    
    """Initiating dust properties calculations based on input parameters"""
    def __init__(self , material , sw , init_dist = 1.0 , size = None , size_range = None):

        self.material = material
        self.sw = sw
        self.init_dist = init_dist
        self.R = init_dist * au

        """Size to beta mapping"""
        if isinstance(size , str):
            self.size = size
            self.r = init_vals[self.size]["r"]
            self.B = init_vals[self.size]["B"][self.material.lower()]
        
        """Entire arrays of sizes and beta values"""
        if size_range is not None:
            self.r , self.B = size_range
        
        self.m0 = size_to_mass(self.r , self.material)

        self.Ytot = self.sputtering_yield()
        self.fsw = self.sw_flux()

        self.V = self.calc_V()
        self.T = self.calc_T()
        self.epsilon = self.eps()

        self.delta = self.V / c
        self.K = self.delta / self.epsilon
        
    """Calculates initial orbital velocity"""   
    def calc_V(self):
        V = np.sqrt((G * m_s * (1 - self.B)) / self.R) #initial orbital velocity, in ms^-1

        return V

    """Calculates initial orbital period"""
    def calc_T(self):
        
        T = np.sqrt(self.R**3 / (G * m_s * (1 - self.B))) #initial orbital period, in s

        return T

    """Calculates solar wind flux using values from Baumann et al., 2020"""
    def sw_flux(self):
        if self.sw == "fast":
            N_sw = 3e6 #particles m^-3, density fast solar wind
            v_sw = 8e5 #ms^-1, velocity fast solar wind
            fsw = N_sw * v_sw #fast solar wind flux, particles/ m^2s^1
    
        elif self.sw == "slow":
            N_sw = 8e6 #particles m^-3, density slow solar wind
            v_sw = 3e5 #ms^-1, velocity slow solar wind
            fsw = N_sw * v_sw #slow solar wind flux, particles/ m^2s^1
    
        elif self.sw == "CME":
            N_sw = 7e7 #particles m^-3, density CME
            v_sw = 5e5 #ms^-1, velocity CME
            fsw = N_sw * v_sw #CME flux, particles/ m^2s^1
    
        return fsw

    """Calculates total sputtering yield"""
    def sputtering_yield(self):
        H = np.sum(sputter[self.material][self.sw]["H"]) #total sputtering yield hydrogen
        He = np.sum(sputter[self.material][self.sw]["He"]) #total sputtering yield helium
        C = np.sum(sputter[self.material][self.sw]["C"]) #total sputtering yield carbon
        O = np.sum(sputter[self.material][self.sw]["O"]) #total sputtering yield oxygen
        N = np.sum(sputter[self.material][self.sw]["N"]) #total sputtering yield nitrogen
        Fe = np.sum(sputter[self.material][self.sw]["Fe"]) #total sputtering yield iron
        Ne = np.sum(sputter[self.material][self.sw]["Ne"]) #total sputtering yield neon
        Mg = np.sum(sputter[self.material][self.sw]["Mg"]) #total sputtering yield magnesium
        Si = np.sum(sputter[self.material][self.sw]["Si"]) #total sputtering yield silicon
        S = np.sum(sputter[self.material][self.sw]["S"]) #total sputtering yield sulfur

        Ytot = H + He + C + O + N + Fe + Ne + Mg + Si + S #summing all total sputtering yields

        return Ytot

    """Calculates sputtering lifetime at 1 AU distance from the Sun, 
    formula from Baumann et al., 2020"""
    def sputtering_lifetime(self):
        if self.material == "silicate":
            rho = rho_s
            mA = mA_S

        elif self.material == "carbon":
            rho = rho_c
            mA = mA_C

        t_sp = (4 * self.r * rho) / (self.fsw * self.Ytot * mA) #sputtering lifetime
         
        self.t_sp = t_sp / yr #sputtering lifetime in units of years

        return self.t_sp

    """Calculates small parameter epsilon_0, mass loss rate"""
    def eps(self):
        if self.material == "silicate":
            rho = rho_s
            mA = mA_S

        elif self.material == "carbon":
            rho = rho_c
            mA = mA_C
        
        eps = (self.fsw * self.Ytot * mA * np.pi * (3 / (4 * np.pi * rho))**(2 / 3) * self.m0**(-1 / 3) 
               * self.T * self.init_dist**(-2))
        
        return eps
    
if __name__ == "__main__":
    # par = dust_properties("silicate" , "slow" , 1 , size = None , size_range = (sil_size , sil_betaval))
    # lower_int = 1 / par.B - np.sqrt((1 - par.B) / par.B**4)
    # upper_int = 1 / par.B + np.sqrt((1 - par.B) / par.B**4)
    # print(lower_int , upper_int)
    par = dust_properties("silicate" , "CME" , 1 , size = "K4")
    print(par.K)
    
    
    