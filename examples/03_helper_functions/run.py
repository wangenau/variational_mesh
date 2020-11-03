#!/usr/bin/env python
'''
Demonstration of helper functions.
'''

from pyscf import dft, gto
from var_mesh import helpers

mol = gto.M(atom='O 0 0 0; H 0 0 0.95691; H 0.95691 0 -0.23987')
mesh = dft.Grids(mol)
mesh.level = 0
mesh.build()

# Plot the mesh with the helper functions
helpers.plot_mesh_3d(mesh=mesh)
helpers.plot_mesh_2d(mesh=mesh, weight=False, plane='xz')
