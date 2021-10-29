import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import sys
sys.path.append("../")
from quelea import *

nx = 217
ny = 133

x0 = 0
x1 = 30 # lambdas
y0 = 0
y1 = 20 # lambdas

xs = np.linspace(x0, x1, nx)
ys = np.linspace(y0, y1, ny)

# 2d array of (x, y, z, t)
coords = np.array( [ [x, y, 0, 0] for x in xs for y in ys ] )
# for map_fields function this should be converted from 2D to 1D array
coords = coords.reshape((4 * nx * ny,))

ftype = 1 # plane wave
a0 = 1 # normalized field amplitude
omega = 1 # frequency
fparam = [a0, 1, 0, 0, 0, 1, 0, 0, omega] # parameters of the plane wave

ex, ey, ez, bx, by, bz = map_fields(coords, ftype, fparam)
# now convert to 2d arrays
ex = ex.reshape((nx, ny))
ey = ey.reshape((nx, ny))
ez = ez.reshape((nx, ny))
bx = bx.reshape((nx, ny))
by = by.reshape((nx, ny))
bz = bz.reshape((nx, ny))
ex = ex.transpose()
ey = ey.transpose()
ez = ez.transpose()
bx = bx.transpose()
by = by.transpose()
bz = bz.transpose()

plt.imshow(ey, cmap = 'RdYlBu', origin = 'lower', extent = [x0, x1, y0, y1])
plt.colorbar()
plt.clim(-a0, a0)

plt.savefig("map_fields.pdf")

