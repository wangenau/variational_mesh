#!/usr/bin/env python
'''
Show the functionality of the mode option.
'''

from pyscf import dft, gto
from var_mesh import var_mesh

mol = gto.M(atom='O 0 0 0; H 0 0 0.95691; H 0.95691 0 -0.23987')
mesh = dft.Grids(mol)

# Call var_mesh with the different mode options
print('mode=PySCF:')
var_mesh(mesh, mode='pyscf')
# If the code does not support different grids for different atom types, the
# precise option has no effect
print('\nmode=ERKALE:')
var_mesh(mesh, precise=True, mode='erkale')
print('\nmode=GAMESS:')
mesh = var_mesh(mesh, precise=False, mode='gamess')
# One can see, that the GAMESS grid differs from the PySCF grid
print('GAMESS grid: %s' % mesh.atom_grid)
