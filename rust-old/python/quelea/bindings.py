# bindings для библиотеки, написанной на Rust

import os
import numpy as np
from ctypes import *

# перекомпилируем библиотеку при импорте этого скрипта
compiled_succ = not os.system("cargo build --release --manifest-path ../../Cargo.toml")

if compiled_succ:

    libfile = "../../target/release/libquelea.so"

    if os.path.isfile(libfile):
        lib = cdll.LoadLibrary(libfile)
        
        #see https://docs.python.org/3/library/ctypes.html
        # dafault lambda value is 1 um
        def drive_particle(s0, t0, dt, n, ftype, fparam, lmbda = 1, rf_scheme = 0):
            nfparam = len(fparam)
            f = lib.drive_particle_binding
            f.restype = POINTER(c_double * (6 * (n + 1)))
            s = (c_double * 6)(*s0)
            fp = (c_double * nfparam)(*fparam)
            v = f(byref(s), c_double(t0), c_double(dt), c_size_t(n), c_int32(ftype), c_size_t(nfparam), byref(fp), c_double(lmbda), c_int(rf_scheme))
            w = np.array(v.contents)
            w = np.reshape(w, (-1, 6))
            return w
        
        def drive_particle_new(s0, t0, dt, n, ftype, fparam):
            nfparam = len(fparam)
            f = lib.drive_particle_new_binding
            f.restype = POINTER(c_double * (6 * (n + 1)))
            s = (c_double * 6)(*s0) # for meaning of * see http://stackoverflow.com/questions/36901/what-does-double-star-and-star-do-for-python-parameters
            fp = (c_double * nfparam)(*fparam)
            v = f(byref(s), c_double(t0), c_double(dt), c_size_t(n), c_int32(ftype), c_size_t(nfparam), byref(fp))
            w = np.array(v.contents)
            w = np.reshape(w, (-1, 6))
            return w

        # Getters for particles
        def x(ps):
            return ps[:,0]
        def y(ps):
            return ps[:,1]
        def z(ps):
            return ps[:,2]
        def ux(ps):
            return ps[:,3]
        def uy(ps):
            return ps[:,4]
        def uz(ps):
            return ps[:,5]
        def g(ps):
            px = ux(ps)
            py = uy(ps)
            pz = uz(ps)
            return np.sqrt(1 + px * px + py * py + pz * pz)
        def t(t0, dt, ps):
            return t0 + np.arange(len(ps)) * dt

        def map_fields(rt, ftype, fparam):
            nrt = len(rt)
            n = nrt // 4
            nfparam = len(fparam)
            f = lib.c_map_fields
            f.restype = POINTER(c_double * (6 * n))
            c_rt = (c_double * nrt)(*rt)
            fp = (c_double * nfparam)(*fparam)
            v = f(c_size_t(n), byref(c_rt), c_int32(ftype), c_size_t(nfparam), byref(fp))
            w = np.array(v.contents)
            ex = w[0::6]
            ey = w[1::6]
            ez = w[2::6]
            bx = w[3::6]
            by = w[4::6]
            bz = w[5::6]
            return (ex, ey, ez, bx, by, bz)

    else:
        print('python/quelea says: file ' + libfile + ' is not found!')

else:
    print('python/quelea says: error while compiling quelea.rs')

