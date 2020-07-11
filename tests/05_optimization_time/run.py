#!/usr/bin/env python
'''
Time the mesh optimization process for different errors.
'''

import matplotlib.pyplot as plt
import numpy as np
import time
from pyscf import dft, gto
from var_mesh import opt_mesh

mol = gto.M(atom='O 0 0 0; H 0 1 0; H 0 0 1')
mesh = dft.Grids(mol)

errors = 10.0**(np.arange(0, -10, -1))
times = []
for i in range(len(errors)):
    start = time.time()
    mesh = opt_mesh(mesh, errors[i])
    end = time.time()
    times.append(end - start)
    print('[%d/%d] %.5Es' % (i + 1, len(errors), end - start))

plt.plot(errors, times)
plt.xlabel('Mesh error')
plt.ylabel('Optimization time [s]')
plt.xscale('log')
plt.show()
