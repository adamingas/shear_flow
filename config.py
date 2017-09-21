import numpy as np
hydro = "yes" #yes or no
noise = "yes" #yes or no
steps = 100
runs = 10
constant = np.linspace(0,5,10).tolist()
time_step = 0.001
yinitial = 0.1
xinitial = 0
zinitial = 0
ra_ratio =  10 #Rmax over radius of sphere
chi = 1 #Thermal energy over max potential energy of FENE spring

