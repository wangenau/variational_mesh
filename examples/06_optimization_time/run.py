#!/usr/bin/env python
'''
Time mesh optimization for different errors.
'''

import matplotlib.pyplot as plt
import numpy as np
import timeit
from pyscf import dft, gto
from var_mesh import var_mesh

mol = gto.M(atom='O 0 0 0; H 0 0 0.95691; H 0.95691 0 -0.23987')
mesh = dft.Grids(mol)

errors = 10.0**(np.arange(-1, -10, -1))
times_false = []
times_true = []
for i in range(len(errors)):
    t1 = timeit.default_timer()
    mesh = var_mesh(mesh, errors[i], precise=False)
    t2 = timeit.default_timer()
    mesh = var_mesh(mesh, errors[i], precise=True)
    t3 = timeit.default_timer()
    times_false.append(t2 - t1)
    times_true.append(t3 - t2)
    print('[%d/%d] Time spent = %f seconds' % (i + 1, len(errors), t3 - t1))

plt.plot(errors, times_false, label='precise=False')
plt.plot(errors, times_true, label='precise=True')
plt.xlabel('Mesh error')
plt.ylabel('Time [s]')
plt.xscale('log')
plt.legend()
plt.show()
