from config import (material_files_bound , au , init_vals , rho_s , rho_c , yr , m_s , sil_size ,
                     sil_betaval , car_betaval_bound , car_size_bound , R , M_ms , M_mc , mA_S , mA_C , size_to_mass)
from sputtering_dict import sputter
from scipy.constants import N_A , G , c
import numpy as np

class dust_properties():
    """Calculating properties of particle given speicified material, solar wind conditions and initial size
    
    Attributes: material (string), "silicate" or "carbon", indicating particle material
                sw (string), "fast", "slow" or "CME", indicating solar wind conditions
                size (string), default: None, "small", "large" or "medium", indicating particle initial size
                size_range (tuple), default: None, else tuple consisting of a set of size values, beta values

                r (float), initial particle size, in m
                B (float), initial beta value, unitless
                m0 (float), initial particle mass, in kg
                Ytot (float), total sputtering yield, unitless
                fsw (float), solar wind flux, in kgm^-3
                V (float), initial orbital velocity, in ms^-1
                T (float), initial orbital period, in s
                epsilon (float), mass loss rate, unitless
                delta (float), PR drag term, unitless
                K (float), mass loss to PR drag ratio, unitless

    Methods: calc_V(), calculates initial orbital velocity
             calc_T(), calculates initial orbital period
             sw_flux(), calculates solar wind flux
             sputtering_yield(), calculates total sputtering yield
             sputtering_lifetime(), calculates sputtering lifetime
             eps(), calculates epsilon
             K_cst_r(), calculates and compares K for all initial r and B values, assuming beta slowly changing

                """
    def __init__(self , material , sw , size = None , size_range = None):
        """Initiating dust properties calculations based on input parameters"""

        self.material = material
        self.sw = sw

        if isinstance(size , str):
            self.size = size
            self.r = init_vals[self.size]["r"]
            self.B = init_vals[self.size]["B"][self.material.lower()]
        
        if size_range is not None:
            self.r , self.B = size_range
        
        self.m0 = size_to_mass(self.r , self.material)

        self.Ytot = self.sputtering_yield()
        self.fsw = self.sw_flux()

        self.V = self.calc_V()
        self.T = self.calc_T()
        self.RAU = R / au #R in AU
        self.epsilon = self.eps()

        self.delta = self.V / c
        self.K = self.delta / self.epsilon
        
        
    def calc_V(self):
        """Calculates initial orbital velocity from given B and R"""
        
        V = np.sqrt((G * m_s * (1 - self.B)) / R) #initial orbital velocity, in ms^-1

        return V

    def calc_T(self):
        """Calculates initial orbital period from given B and R"""
        T = np.sqrt(R**3 / (G * m_s * (1 - self.B))) #initial orbital period, in s

        return T

    def sw_flux(self):
        """Calculates solar wind flux based on user input solar wind conditions"""
    
        if self.sw == "fast":
            N_sw = 3e6 #particles m^-3, density fast solar wind
            v_sw = 8e5 #ms^-1, velocity fast solar wind
            fsw = N_sw * v_sw #solar wind flux, particles m^-2s^-1
    
        elif self.sw == "slow":
            N_sw = 8e6 #particles m^-3, density slow solar wind
            v_sw = 3e5 #ms^-1, velocity slow solar wind
            # v_sw = 8e5
            fsw = N_sw * v_sw #solar wind flux, particles m^-2s^-1
    
        elif self.sw == "CME":
            N_sw = 7e7 #particles m^-3, density slow solar wind
            v_sw = 5e5 #ms^-1, velocity slow solar wind
            fsw = N_sw * v_sw #solar wind flux, particles m^-2s^-1
    
        return fsw

    #Total sputtering yield
    def sputtering_yield(self):
        """Calculates total sputtering yield based on input material and solar wind conditions"""
    
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

        Ytot = H + He + C + O + N + Fe + Ne + Mg + Si + S #total sputtering yield 

        return Ytot

    def sputtering_lifetime(self):
        """Calculates sputtering lifetime from input material and solar wind conditions, 1 AU distance from the Sun"""

        if self.material == "silicate":
            M_m = M_ms
            rho = rho_s
            mA = mA_S

        elif self.material == "carbon":
            M_m = M_mc
            rho = rho_c
            mA = mA_C

        t_sp = (4 * self.r * rho) / (self.fsw * self.Ytot * mA) #sputtering lifetime
         
        self.t_sp = t_sp / yr

        return self.t_sp

    """calculates small parameter epsilon, using constants from mass calcs"""
    def eps(self):
        """input: m (float), default M, mass of particle in kg
                sw (string), default: slow, options fast and CME
                species (string), default: all, else one of the elements H, He, C, O, N, Fe, Ne, Mg, Si, S

        returns: eps (float), epsilon parameter"""
    
        if self.material == "silicate":
            M_m = M_ms
            rho = rho_s
            mA = mA_S

        elif self.material == "carbon":
            M_m = M_mc
            rho = rho_c
            mA = mA_C
        
        eps = self.fsw * self.Ytot * mA * np.pi * (3 / (4 * np.pi * rho))**(2 / 3) * self.m0**(-1 / 3) * self.T * self.RAU**(-2)
        
        return eps

    def K_cst_r(self):
        Kcst = (1 - self.B)**2 / 6

        tol = 1e-5
        
        vals = np.zeros((3 , len(Kcst)))

        for i in range(len(self.K)):
            
            if self.K[i] - Kcst[i] < tol:
                vals[: , i] = [self.K[i] , self.B[i] , self.r[i]]

        return vals
    
if __name__ == "__main__":
    
    par = dust_properties("silicate" , "slow" , "large")
    # par = dust_properties("silicate" , "CME" , size = "large")
    
    print(par.K)
    
    
    