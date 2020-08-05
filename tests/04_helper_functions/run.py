#!/usr/bin/env python
'''
Demonstration of helper functions.
'''

from pyscf import dft, gto
from var_mesh import helpers

mol = gto.M(atom='O 0 0 0; H 0 1 0; H 0 0 1')
mesh = dft.Grids(mol)
mesh.level = 0
mesh.build()

helpers.plot_mesh_3d(mesh=mesh)
helpers.plot_mesh_2d(mesh=mesh, weight=False, plane='yz')
