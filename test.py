import numpy as np
import matplotlib.pyplot as plt

rk = np.load("Files/rk45_t5_large_silicate_slowsw_realbeta_test.npz")
x1 , y1 , vx1 , vy1 , m1 , b1 = [rk[k] for k in ("x" , "y" , "vx" , "vy" , "m" , "b")]

plt.plot(x1 , y1)
plt.show()