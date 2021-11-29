#see Fig. 6(a) from [Samsonov 2018, PRA]

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import sys
sys.path.append("../")
from quelea import *

a0 = 1500
wavelength = 1
width = 2.5

gamma0 = 300
#initial momenta and coords of the electron
s0 = [0, 0, 0, 0.8 * gamma0, 0, 0.6 * gamma0]

t_end = 50 # lambda / c
dt = 0.001
nt = int( np.round(t_end / dt) )

traj = drive_particle(s0, 0, dt, nt, 15, [width, a0], wavelength, 1) #(r: (f64,f64, f64), r0: f64, t: f64, a0: f64)
plt.plot(x(traj), g(traj))
plt.xlabel(r"$x$")
plt.ylabel(r"$g$")
#plt.plot(t(0, dt, traj), g(traj))

plt.savefig("gauss.pdf")
