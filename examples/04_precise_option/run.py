#!/usr/bin/env python
'''
Show the impact of the precise error option.
'''

from pyscf import dft, gto
from timeit import default_timer
from var_mesh import var_mesh

# Set an error threshold
error = 1e-8
mol = gto.M(atom='O 0 0 0; H 0 0 0.95691; H 0.95691 0 -0.23987')
mf = dft.RKS(mol)

# Time mesh generation with the precise option disabled
print('precise=False:')
start = default_timer()
mesh = var_mesh(mf.grids, thres=error, precise=False)
end = default_timer()
print('Time spent = %f seconds' % (end - start))
print('Mesh points = %d\n' % len(mesh.coords))

# Time mesh generation with the precise option enabled
print('precise=True:')
start = default_timer()
mesh = var_mesh(mf.grids, thres=error, precise=True)
end = default_timer()
print('Time spent = %f seconds' % (end - start))
print('Mesh points = %d' % len(mesh.coords))
