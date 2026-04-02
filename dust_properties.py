from config import sil_beta , car_beta , init_vals , rho , yr , m_s , u , m_Mg , m_Si , m_O , m_Fe , m_C , R , M_ms , M_mc , mA_S , mA_C
from AstronomicalSilicate_modified import sputter
from scipy.constants import N_A , G , c
import numpy as np

class dust_properties():

    def __init__(self , material , sw , species , size):

        self.material = material
        self.sw = sw
        self.species = species #solar wind elements sputtering
        self.size = size
        self.r = init_vals[self.size]["r"]
        self.B = init_vals[self.size]["B"][self.material.lower()]
        self.m0 = self.m()

        self.V = np.sqrt((G * m_s * (1 - self.B)) / R) #initial angular velocity, scaled formula
        self.T = round(np.sqrt(R**3 / (G * m_s * (1 - self.B))))

        self.Ytot = self.sputtering_yield()
        self.fsw = self.sw_flux()

        self.delta = self.V / c

    def m(self):
        mass = 4 / 3 * np.pi * rho * (self.r)**3

        return mass

    #Solar wind flux
    def sw_flux(self):
        """input: sw (string), solar wind conditions, default slow, options fast, slow, CME
        
            returns: fsw (array), solar wind flux"""
    
        if self.sw == "fast":
            N_sw = 3e6 #m^-3 density fast solar wind
            v_sw = 8e5 #ms^-1, velocity fast solar wind
            fsw = N_sw * v_sw #solar wind flux, m^-2s^-1
    
        elif self.sw == "slow":
            N_sw = 8e6 #m^-3, density slow solar wind
            v_sw = 3e5 #ms^-1, velocity slow solar wind
            fsw = N_sw * v_sw #solar wind flux, m^-2s^-1
    
        elif self.sw == "CME":
            N_sw = 7e7 #m^-3, density slow solar wind
            v_sw = 5e5 #ms^-1, velocity slow solar wind
            fsw = N_sw * v_sw #solar wind flux, m^-2s^-1
    
        return fsw

    #Total sputtering yield
    def sputtering_yield(self):
        """input: sw (string), solar wind conditions, default slow, options fast, slow, CME
                species (string), solar wind elements, default all, else one of H, He, C, O, N, Fe, Ne, Mg, Si, S
        
            returns: Ytot (array), total sputtering yield"""
    
        if self.species == "all":
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
    
        else:
            spec = np.sum(sputter[self.material][self.sw][self.species]) #total sputtering yield specified species

            Ytot = spec #total sputtering yield

        return Ytot

    """sputtering lifetime calcs"""
    def sputtering_lifetime(self):
        """input: r0 (float), initial dust radius
                fsw (float), solar wind flux
                Ytot (float), total sputtering yield
                M_m (float), molar mass of the material

        returns: t_sp (float), sputtering lifetime"""

        if self.material == "silicate":
            M_m = M_ms

        elif self.material == "carbon":
            M_m = M_mc

        t_sp = (4 * self.r * N_A * rho) / (self.fsw * self.Ytot * M_m) #sputtering lifetime
        self.t_sp = t_sp / yr

        return self.t_sp

    """calculates small parameter epsilon, using constants from mass calcs"""
    def eps(self):
        """input: m (float), default M, mass of particle in kg
                sw (string), default: slow, options fast and CME
                species (string), default: all, else one of the elements H, He, C, O, N, Fe, Ne, Mg, Si, S

        returns: eps (float), epsilon parameter"""
    
        if self.material == "silicate":
            mA = mA_S
    
        elif self.material == "carbon": 
            mA = mA_C
    
        eps = self.fsw * self.Ytot * mA * np.pi * (3 / (4 * np.pi * rho))**(2 / 3) * self.m0**(-1 / 3) * self.T 

        return eps