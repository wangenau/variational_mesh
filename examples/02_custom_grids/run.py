#!/usr/bin/env python
'''
Use custom grids for mesh optimization.
'''

from pyscf import dft, gto
from var_mesh import gen_mesh

# Set custom radial grids
# Setting custom grids only for a subset of atomic species also works
gen_mesh.rad = {'H': list(range(10, 110, 10)),
                'O': list(range(20, 220, 20))}
# Set custom angular grids
# These grids have to follow the Lebedev order
gen_mesh.ang = {'H': gen_mesh.ang_grids[15:20],
                'O': gen_mesh.ang_grids[20:25]}

mol = gto.M(atom='O 0 0 0; H 0 0 0.95691; H 0.95691 0 -0.23987')
mesh = dft.Grids(mol)
# Get the maximal output with a larger verbosity level
mesh.verbose = 5
# Change mesh parameters before calling var_mesh
mesh.prune = None
mesh = gen_mesh.var_mesh(mesh, thres=1e-7)
print('Mesh points = %d' % len(mesh.coords))
