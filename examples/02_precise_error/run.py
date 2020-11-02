#!/usr/bin/env python
'''
Impact of the precise error option.
'''

import timeit
from pyscf import dft, gto
from var_mesh import var_mesh

thres = 1e-8
mol = gto.M(atom='O 0 0 0; H 0 0 0.95691; H 0.95691 0 -0.23987')
mf = dft.RKS(mol)

start = timeit.default_timer()
var_mesh(mf.grids, thres, precise=False)
end = timeit.default_timer()
print('\'precise=False\' = %f seconds' % (end - start))

start = timeit.default_timer()
var_mesh(mf.grids, thres, precise=True)
end = timeit.default_timer()
print('\'precise=True\'  = %f seconds' % (end - start))
