#!/usr/bin/env python
'''
Minimal example to use the variational mesh.
'''

from pyscf import dft, gto
from var_mesh import var_mesh

# Create a mol object for H2O
mol = gto.M(atom='O 0 0 0; H 0 0 0.95691; H 0.95691 0 -0.23987')

# Create an object for a restricted Kohn-Sham SCF calculation
mf = dft.RKS(mol)

# Overwrite the grids object by calling the var_mesh function
mf.grids = var_mesh(mf.grids)

# Run the calculation
mf.kernel()
