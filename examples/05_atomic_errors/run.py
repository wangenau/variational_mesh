#!/usr/bin/env python
'''
Generate mesh errors for an atom with standard PySCF grid levels.
'''

import matplotlib.pyplot as plt
from pyscf import dft, gto
from var_mesh.gen_mesh import mesh_error

mol = gto.M(atom='H 0 0 0', spin=1, verbose=0)
mf = dft.UKS(mol)
mf.max_cycle = 0

errors = []
for i in range(10):
    print('Grid level = %d' % i)
    mf.grids.level = i
    mf.grids.build()
    mf.kernel()
    error = mesh_error(mf)
    errors.append(error)
    print('Mesh error = %.5e' % error)

plt.plot(range(10), errors)
plt.xlabel('Grid level')
plt.ylabel('Mesh error')
plt.yscale('log')
plt.show()
