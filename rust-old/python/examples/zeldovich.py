#see Fig. 6(a) from [Samsonov 2018, PRA]

# Numerical solution for the Zeldovich problem in E_x,y rotating around the Z axis. Initially
# electric field is opposite to the X axis.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.constants as cns


import sys
sys.path.append("../")
from quelea import *

a0 = 400
wavelength = 1 # um
fparams = [a0, 0, 0, 1, -1, 0, 0]

#initial momenta and coords of the electron
s0 = [0, 0, 0, 0, 0, 0]

t_end = 30 # lambda / c, для `phase_space` предполагается, что это - целое число.
dt = 0.01

n_traj = 4000

# фазовое пространство в "целые" моменты времени `nT`, когда вектор E направлен вдоль оси x.
def phase_space(n_traj, a0, wavelength, s0, nT, dt):
    pxs = np.empty(n_traj)
    pys = np.empty(n_traj)
    A = 2/3 * cns.alpha * cns.hbar * ( cns.c * 2 * cns.pi / (wavelength * 1e-6) ) / (cns.m_e * cns.c**2) * a0**3
    p = np.sqrt( ( np.sqrt( 4 * A**2 + 1 ) - 1 ) / ( 2 * A**2 ) )
    psi = np.arcsin(p)
    
    for i in range(n_traj):
        traj = drive_particle( s0, 0, dt, int( int(nT) * 2 * np.pi / dt ), 3, fparams, wavelength, 1)
        pxs[i] = ux(traj)[-1] / a0
        pys[i] = uy(traj)[-1] / a0

    pxmin = np.min([0, np.min(pxs)])
    pymax = np.max([0, np.max(pys)])
    plt.hist2d(pxs, pys, bins = [50, 50], range = [[pxmin, np.max(pxs)], [np.min(pys), pymax]], cmap = 'viridis')
    plt.plot([0, p * np.cos(psi)], [0, -p * np.sin(psi)], 'w:')
    plt.plot(p * np.cos(psi), -p * np.sin(psi), 'wo')
    plt.plot(0, 0, 'w*')
    ax = plt.gca()
    ax.set_aspect(1)
    plt.xlabel(r"$p_x$")
    plt.ylabel(r"$p_y$")
    plt.savefig("zeldovich_phase_space_" + str(a0) + "_" + str(wavelength) + ".png")

# trajectrory of a single electron in px, py space
def traj(t_end):
    traj = drive_particle(s0, 0, dt, int( np.round(t_end * 2 * np.pi / dt) ), 3, fparams, wavelength, 1)
    px = ux(traj) / a0
    py = uy(traj) / a0
    plt.plot(px, py, 'b-')
    plt.savefig("zeldovich_traj_" + str(a0) + "_" + str(wavelength) + ".png")

#    u_x = -np.power(p,3)*A*np.cos(phi) - np.power(p,2)*np.sin(phi)
#    u_y = -np.power(p,3)*A*np.sin(phi) + np.power(p,2)*np.cos(phi)
#
#    #plt.hist2d(abscissa, ordinate, bins=[100,100], cmap=plt.cm.jet, label=r"a_0 = 800, \lambda = 0.01")
#    plt.plot(u_x, u_y)
#    plt.xlabel(r"$p_x$")
#    plt.ylabel(r"$p_y$")
#    #plt.plot(t(0, dt, traj), g(traj))
#    plt.savefig("zeldovich_" + str(a0) + "_" + str(wavelength) + ".png")

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
    phase_space(n_traj, a0, wavelength, s0, t_end, dt)
    #traj(t_end)
    
if __name__ == "__main__":
    main()
