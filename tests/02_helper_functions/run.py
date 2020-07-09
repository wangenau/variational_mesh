#!/usr/bin/env python
'''
Usage of helper functions.
'''

from pyscf import dft, gto
from var_mesh.helpers import *

mol = gto.M(atom='H -1.08 0 0; C 0 0 0; N 1.15 0 0')
mesh = dft.Grids(mol)
mesh.level = 0
mesh.build()

plot_mesh_3d(mesh=mesh)
plot_mesh_2d(mesh=mesh, weight=False, plane='xz')
plot_combs(mol=mol, level=5)
