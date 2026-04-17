import numpy as np
import scipy.integrate as it
from dust_properties import dust_properties
import matplotlib.pyplot as plt

def C0(B , K , beta , time , epsilon):
    cst = - 4 * B * K / (1 - B)**3
    var = beta * (1 - B * beta)**4
    b_int = it.cumulative_trapezoid(var , time , initial = 0)
    terms = 1.0 + cst * b_int * epsilon
    # invalid = np.where(terms <= 0)
    # terms[invalid] = 0.0000001

    C0_fourth = terms**(1 / 4)

    return C0_fourth

if __name__ == "__main__":
    par = dust_properties("silicate" , "slow" , "all" , "large")
    res = np.load("Files/rk45_t6_large_silicate_slowsw.npz")
    x , y , vx , vy , m , bnum , t = [res[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b", "t")]
    c0func = C0(par.B , par.K , bnum , t , par.epsilon)
    print(c0func)
    # plt.plot(t , c0func)
    # plt.show()
    

    