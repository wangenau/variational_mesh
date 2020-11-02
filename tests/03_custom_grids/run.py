#!/usr/bin/env python
'''
Use custom grids for the mesh optimization.
'''

from pyscf import dft, gto
from var_mesh import gen_mesh

gen_mesh.rad = {'H': list(range(10, 110, 10)),
                'O': list(range(20, 220, 20))}
gen_mesh.ang = {'H': gen_mesh.ang_grids[15:20],
                'O': gen_mesh.ang_grids[20:25]}

mol = gto.M(atom='O 0 0 0; H 0 0 0.95691; H 0.95691 0 -0.23987')
mol.verbose = 5
mesh = dft.Grids(mol)
mesh = gen_mesh.var_mesh(mesh, 1e-7)
print(mesh.coords.shape)
