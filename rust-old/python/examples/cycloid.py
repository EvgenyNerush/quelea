#see Fig. 6(a) from [Samsonov 2018, PRA]

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import sys
sys.path.append("../")
from quelea import *

a0 = 50.3
wavelength = 1.1e-3 # 1.1 nm

gamma0 = 63
#initial momenta and coords of the electron
s0 = [0, 0, 0, 0.8 * gamma0, 0, 0.6 * gamma0]

t_end = 21.47 # lambda / c
dt = 0.02
nt = int( np.round(t_end / dt) )

traj = drive_particle(s0, 0, dt, nt, 0, [a0, 0], wavelength, 0)

plt.plot(x(traj), y(traj))
#plt.plot(t(0, dt, traj), g(traj))

plt.savefig("cycloid.pdf")
