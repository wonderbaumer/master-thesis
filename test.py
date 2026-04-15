import numpy as np
import matplotlib.pyplot as plt

"""plots B values for masses corresponding to a range of initial particle sizes, assuming beta=Frad/G, not real beta curve
and epsilon calculated same range of masses"""
def eps_init_beta():
    """input: none
       
       returns: none
       
       tror ikke denne egt funker som den skal, regner bare V en gang, ikke for hver B"""
    
    material = ["silicate" , "carbon"]
                
    
    ls = {"silicate" : "-" , 
                "carbon" : "--"}
    
    sw_conds = ["slow" , "fast" , "CME"]
    eps_vals = []

    for mat in material:
        for sw in sw_conds:
            par = dust_properties(mat , sw , "all" , None , b_func = True)
            epsilon , newB = par.update_bfunc()
            eps_vals.append({"material" : mat ,
                             "sw_cond" : sw ,
                             "eps" : epsilon ,
                             "betaval" : newB})
    
    for i in eps_vals:
        label = f"{i['material']} for {i['sw_cond']}"
        vals = i["eps"]
        b_init_vals = i["betaval"]
        plt.plot(b_init_vals , vals , linestyle = ls[i['material']] , label = label)

    plt.xlabel(r"$B$")
    plt.ylabel(r"${\epsilon}$")
    plt.yscale("log")
    plt.title(r"${\epsilon}$ vs ${B}$, corresponding to size range $1~\mathrm{nm} \text{–} 50~\mu\mathrm{m}$, silicate and carbon")
    plt.legend(loc = "lower right")
    plt.savefig("Plots/eps_vs_beta.png" , dpi = 300 , bbox_inches = 'tight')
    plt.show()

       """
    def du(self , t , u):
        #c0 = c0[0]
        #c0_safe = max(c0, 1e-12)
        beta = self.beta_func(t)
        tot = np.array([-4 * self.B * self.K * beta * (1 - beta * self.B)**4 / ((1 - self.B)**3)])

        return tot
    
    def solve_c0(self):

        ttot = (self.time[0] , self.time[-1])

        sol = solve_ivp(self.du , ttot , [self.C0_0**4] , method = "Radau" , t_eval = self.time , rtol=1e-9 , atol=1e-12)

        self.u = sol.y[0]
        self.beta_vals = self.beta_func(self.time)
        self.c0 = self.u**(1 / 4)

        return self.c0
    """

    #self.beta_func = pchip(self.time , self.barr , extrapolate = False)
        # self.C0_0 = 1.0
        # self.c0 = None
        # self.c0_func = None
        # self.beta_vals = None
        #print(self.epsilon / 2.6031255888589655) #1.9279315785095778e-06)
    
    """
    def betahat_analytical(self):
        bvals = self.epsilon * self.time / 2.6031255888589655 + 0.9996231879878003
        return bvals
    
    """ 