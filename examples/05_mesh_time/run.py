#!/usr/bin/env python
'''
Time the mesh generation for different error thresholds.
'''

import matplotlib.pyplot as plt
import numpy as np
from pyscf import dft, gto
from timeit import default_timer
from var_mesh import var_mesh

mol = gto.M(atom='O 0 0 0; H 0 0 0.95691; H 0.95691 0 -0.23987')
mesh = dft.Grids(mol)
mesh.verbose = 0

# Set up different mesh errors and generate a mesh for it
errors = 10.0**(np.arange(-1, -9, -1))
times_false = []
times_true = []
for i in range(len(errors)):
    print('[%d/%d]' % (i + 1, len(errors)))
    print('Error threshold = %.0e' % errors[i])
    t1 = default_timer()
    mesh = var_mesh(mesh, thres=errors[i], precise=False)
    t2 = default_timer()
    mesh = var_mesh(mesh, thres=errors[i], precise=True)
    t3 = default_timer()
    times_false.append(t2 - t1)
    times_true.append(t3 - t2)
    print('Time spent = %f seconds' % (t3 - t1))

# Plot everything
plt.plot(errors, times_false, label='precise=False')
plt.plot(errors, times_true, label='precise=True')
plt.xlabel('Mesh error')
plt.ylabel('Time [s]')
plt.xscale('log')
plt.gca().invert_xaxis()
plt.legend()
plt.show()
