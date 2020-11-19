#!/usr/bin/env python3
'''
Show usage with the PyFLOSIC package. Requires ase and PyFLOSIC v1.0.X.
Python3 is required for PyFLOSIC. This example is based on:
https://github.com/pyflosic/pyflosic/tree/master/examples/basic_calculations
'''

try:
    from ase.io import read
except ImportError:
    raise SystemExit('The package \'ase\' is required for this example!')
try:
    from flosic_os import ase2pyscf, xyz_to_nuclei_fod
    from flosic_scf import FLOSIC
except ImportError:
    raise SystemExit('The package \'pyflosic\' is required for this example!')
from pyscf import gto
from var_mesh import var_mesh

# Set up calculation details
molecule = read('H2.xyz')
geo, nuclei, fod1, fod2, included = xyz_to_nuclei_fod(molecule)
mol = gto.M(atom=ase2pyscf(nuclei), basis='6-311++Gss', spin=0, charge=0)
sic_object = FLOSIC(mol, xc='lda,pw', fod1=fod1, fod2=fod2, ham_sic='HOO')
sic_object.max_cycle = 300
sic_object.conv_tol = 1e-7

# The default example uses a grid level of 4, save its size
sic_object.grid_level = 4
mesh_size = len(sic_object.calc_uks.grids.coords)
sic_object.calc_uks.grids = var_mesh(sic_object.calc_uks.grids)

# Start the calculation
total_energy_sic = sic_object.kernel()
homo_flosic = sic_object.homo_flosic

# Display mesh sizes and energy values
print('Mesh size before: %d' % mesh_size)
print('Mesh size after: %d' % len(sic_object.calc_uks.grids.coords))
print('Total energy of H2 (FLO-SIC SCF): %0.5f (should be %0.5f)' %
      (total_energy_sic, -1.18118689724))
print('HOMO energy eigenvalue of H2 (FLO-SIC SCF): %0.5f (should be %0.5f)' %
      (homo_flosic, -0.623425516328))
