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