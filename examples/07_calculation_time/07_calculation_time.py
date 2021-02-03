#!/usr/bin/env python
'''
Time the mesh generation with the respective SCF calculation.
'''

import matplotlib.pyplot as plt
import numpy as np
from pyscf import dft, gto
from timeit import default_timer
from var_mesh import var_mesh

mol = gto.M(atom='O 0 0 0; H 0 0 0.95691; H 0.95691 0 -0.23987')
mf = dft.RKS(mol)
mf.verbose = 0
mf.grids.verbose = 0

# Set up different mesh errors and generate a mesh for it
errors = 10.0**(np.arange(-1, -9, -1))
time_mesh = []
time_scf = []
for i in range(len(errors)):
    print('[%d/%d]' % (i + 1, len(errors)))
    print('Error threshold = %.0e' % errors[i])
    t1 = default_timer()
    mf.grids = var_mesh(mf.grids, thres=errors[i], precise=True)
    t2 = default_timer()
    mf.kernel()
    t3 = default_timer()
    time_mesh.append(t2 - t1)
    time_scf.append(t3 - t2)
    print('Time spent = %f seconds' % (t3 - t1))

# Plot everything
plt.plot(errors, time_mesh, label='VarMesh')
plt.plot(errors, time_scf, label='SCF')
plt.xlabel('Mesh error')
plt.ylabel('Time [s]')
plt.xscale('log')
plt.gca().invert_xaxis()
plt.legend()
plt.show()
