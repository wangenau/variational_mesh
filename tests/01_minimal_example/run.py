#!/usr/bin/env python
'''
Minimal example to generate a variational mesh.
'''

from pyscf import dft, gto
from var_mesh import gen_mesh

mol = gto.M(atom='O 0 0 0; H 0 1 0; H 0 0 1')
mf = dft.RKS(mol)
mf.grids = gen_mesh(mf.grids, 1e-5)
mf.kernel()
