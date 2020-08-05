#!/usr/bin/env python
'''
Generate total energies error and mesh errors.
'''

import matplotlib.pyplot as plt
import numpy as np
from pyscf import dft, gto
from var_mesh import opt_mesh

mol = gto.M(atom='O 0 0 0; H 0 1 0; H 0 0 1')
mf = dft.RKS(mol)

errors = 10.0**(np.arange(-6, -11, -1))
etot = []
for i in range(len(errors)):
    print('[%d/%d]' % (i + 1, len(errors)))
    mf.grids = opt_mesh(mf.grids, errors[i])
    mf.kernel()
    etot.append(mf.e_tot)

fig, ax1 = plt.subplots()
ax1.set_xlabel('Steps')
ax1.set_ylabel('Mesh error', color='blue')
ax1.plot(range(len(errors)), errors, color='blue')
ax2 = ax1.twinx()
ax2.set_ylabel('Total energy [Ha]', color='orange')
ax2.plot(range(len(etot)), etot, color='orange')
plt.tight_layout()
plt.show()
