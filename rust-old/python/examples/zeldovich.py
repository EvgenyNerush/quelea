#see Fig. 6(a) from [Samsonov 2018, PRA]

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.constants as cns


import sys
sys.path.append("../")
from quelea import *

a0 = 1500
wavelength = 1

#initial momenta and coords of the electron
s0 = [0, 0, 0, 0, 0, 0]

t_end = 100 # lambda / c
dt = 0.01
nt = int( np.round(t_end / dt) )

traj_n = 1000

def phase_portrate(traj_n, a0, wavelength, s0, nt, dt):
    abscissa = np.empty(traj_n)
    ordinate = np.empty(traj_n)
    A = 2/3 *cns.alpha*cns.hbar/np.power(cns.c,2)*a0**3
    p = (np.sqrt(np.sqrt(4*np.power(A,2)+1)-1))/(2*np.power(A,2))
    
    for i in range(traj_n):
        traj = drive_particle(s0, 0, dt, nt, 3, [a0, 0, 0, 1, a0, 0, 0], wavelength, 1) 
        abscissa[i] = ux(traj)[-1]
        ordinate[i] = uy(traj)[-1]

    traj = drive_particle(s0, 0, dt, nt, 3, [a0, 0, 0, 1, a0, 0, 0], wavelength, 1)
    phi = t(0, dt, traj)

    u_x = -np.power(p,3)*A*np.cos(phi) - np.power(p,2)*np.sin(phi)
    u_y = -np.power(p,3)*A*np.sin(phi) + np.power(p,2)*np.cos(phi)

    #plt.hist2d(abscissa, ordinate, bins=[100,100], cmap=plt.cm.jet, label=r"a_0 = 800, \lambda = 0.01")
    plt.plot(u_x, u_y)
    plt.xlabel(r"$p_x$")
    plt.ylabel(r"$p_y$")
    #plt.plot(t(0, dt, traj), g(traj))
    plt.savefig("zeldovich_" + str(a0) + "_" + str(wavelength) + ".png")

def spiral_plot(a0, wavelength, s0, nt, dt):
    traj = drive_particle(s0, 0, dt, nt, 3, [a0, 0, 0, 1, a0, 0, 0], wavelength, 0) #LL force only
    phi = t(0, dt, traj) #angle between  abscissa and E-parallel axis
    u_parallel = ux(traj)*np.cos(phi) + uy(traj)*np.sin(phi)
    u_perp = uy(traj)*np.cos(phi) - ux(traj)*np.sin(phi)

    plt.plot(u_parallel, u_perp)
    plt.xlabel(r'$p_\parallel$')
    plt.ylabel(r'$p_\perp$')
    plt.savefig("pp_ort_" + str(a0) + "_" + str(wavelength) + ".png")



def main():
    #for a0 in []:
    #for wavelength in []:
    #spiral_plot(a0, wavelength, s0, nt, dt)
    phase_portrate(traj_n, a0, wavelength, s0, nt, dt)
    
if __name__ == "__main__":
    main()