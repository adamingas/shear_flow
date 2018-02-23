import numpy as np

noise = "yes" #yes or no
steps = 100000
runs = 2
constant = np.linspace(0,10,21).tolist()
time_step = 0.01
yinitial = 0.1
xinitial = 0
zinitial = 0
ar_ratio =  0.1 #Rmax over radius of sphere
chi = [10**(-3)] #Thermal energy over max potential energy of FENE spring
if ar_ratio == 0:
	hydro = "no"
else:
	hydro = "yes" 
