#!/usr/bin/env python
'''
Generate mesh errors for an atom with standard PySCF grid levels.
'''

from pyscf import dft, gto
from var_mesh.gen_mesh import mesh_error
import matplotlib.pyplot as plt

mol = gto.M(atom='C 0 0 0', verbose=0)
mf = dft.RKS(mol)
mf.max_cycle = 0

errors = []
for i in range(10):
    print('Grid level: %d' % i)
    mf.grids.level = i
    mf.grids.build()
    mf.kernel()
    error = mesh_error(mf)
    errors.append(error)
    print('Mesh error: %.5E' % error)

plt.plot(range(10), errors)
plt.yscale('log')
plt.xlabel('Grid level')
plt.ylabel('Mesh error')
plt.show()
