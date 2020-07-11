#!/usr/bin/env python
'''
Minimal example to use the variational mesh.
'''

from pyscf import dft, gto
from var_mesh import opt_mesh

mol = gto.M(atom='O 0 0 0; H 0 1 0; H 0 0 1')
mf = dft.RKS(mol)
mf.grids = opt_mesh(mf.grids, 1e-5)
mf.kernel()