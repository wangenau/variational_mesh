#!/usr/bin/env python
'''
Time mesh optimization and SCF calculation.
'''

import matplotlib.pyplot as plt
import numpy as np
import timeit
from pyscf import dft, gto
from var_mesh import var_mesh

mol = gto.M(atom='O 0 0 0; H 0 0 0.95691; H 0.95691 0 -0.23987')
mf = dft.RKS(mol)

errors = 10.0**(np.arange(-1, -10, -1))
time_mesh = []
time_scf = []
for i in range(len(errors)):
    t1 = timeit.default_timer()
    mf.grids = var_mesh(mf.grids, errors[i], precise=True)
    t2 = timeit.default_timer()
    mf.kernel()
    t3 = timeit.default_timer()
    time_mesh.append(t2 - t1)
    time_scf.append(t3 - t2)
    print('[%d/%d] Time spent = %f seconds' % (i + 1, len(errors), t3 - t1))

plt.plot(range(len(time_mesh)), time_mesh, label='Mesh')
plt.plot(range(len(time_scf)), time_scf, label='SCF')
plt.xlabel('Steps')
plt.ylabel('Time [s]')
plt.legend()
plt.show()
