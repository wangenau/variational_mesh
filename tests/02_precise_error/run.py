#!/usr/bin/env python
'''
Test the impact of the precise error option.
'''

import timeit
from pyscf import dft, gto
from var_mesh import opt_mesh

mol = gto.M(atom='O 0 0 0; H 0 1 0; H 0 0 1')
mf = dft.RKS(mol)

start = timeit.default_timer()
mf.grids = opt_mesh(mf.grids, 1e-11, precise=False)
end = timeit.default_timer()
print('precise=False: %.5Es' % (end - start))

start = timeit.default_timer()
mf.grids = opt_mesh(mf.grids, 1e-11, precise=True)
end = timeit.default_timer()
print('precise=True:  %.5Es' % (end - start))
