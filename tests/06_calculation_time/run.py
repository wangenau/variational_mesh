#!/usr/bin/env python
'''
Time the mesh optimization and SCF calculation.
'''

import matplotlib.pyplot as plt
import numpy as np
import time
from pyscf import dft, gto
from var_mesh import opt_mesh

mol = gto.M(atom='O 0 0 0; H 0 1 0; H 0 0 1')
mf = dft.RKS(mol)

errors = 10.0**(np.arange(0, -10, -1))
time_mesh = []
time_etot = []
for i in range(len(errors)):
    t1 = time.time()
    mf.grids = opt_mesh(mf.grids, errors[i])
    t2 = time.time()
    mf.kernel()
    t3 = time.time()
    time_mesh.append(t2 - t1)
    time_etot.append(t3 - t2)
    print('[%d/%d] %.5Es' % (i + 1, len(errors), t3 - t1))

plt.plot(range(len(time_mesh)), time_mesh, label='Mesh')
plt.plot(range(len(time_etot)), time_etot, label='Energy')
plt.xlabel('Steps')
plt.ylabel('Time [s]')
plt.legend()
plt.show()
