import matplotlib.pyplot as plt
import numpy as np
from pyscf import dft, gto
from pyscf.dft import radi


def gen_partition_equidist(mol, atom_grids_tab, radii_adjust=None,
                           atomic_radii=radi.BRAGG_RADII, becke_scheme=None):
    '''Copy of equidistant_mesh.gen_partition_equidist but nsteps is global.'''
    atm_coords = np.asarray(mol.atom_coords())
    coords_all = []
    weights_all = []
    for ia in range(mol.natm):
        chg = gto.charge(mol.atom_symbol(ia))
        rad = atomic_radii[chg]
        steps = np.linspace(-rad, rad, nsteps)
        coords = []
        for x in steps:
            for y in steps:
                for z in steps:
                    coords.append([x, y, z])
        coords = np.vstack(coords)
        coords += atm_coords[ia]
        vol = (2 * rad / nsteps)**3
        weights = np.full(len(coords), vol)
        coords_all.append(coords)
        weights_all.append(weights)
    return np.vstack(coords_all), np.hstack(weights_all)


mol = gto.M(atom='O 0, 0, 0; H 0, 1, 0; H 0, 0, 1', basis='sto-3g')
mf = dft.RKS(mol).set(xc='lda,pw')

# calculate the total energy for a rising number of mesh points
e_tot = []
grids = []
mf.grids.gen_partition = gen_partition_equidist
for i in range(10):
    print('Grid level: {}'.format(i))
    # necessary for mf.kernel to start a new calculation
    mf.grids.level = i
    nsteps = (i + 1) * 10
    mf.kernel()
    e_tot.append(mf.e_tot)
    grids.append(nsteps**3)

# plot energies per total grid points
fig, ax = plt.subplots()
plt.xlabel('grid points')
plt.ylabel('energies [Ha]')
plt.plot(grids, e_tot, label='equidistant mesh')
plt.xscale('log')
# pyscf energy value for H2O, calculated at grid level 9
ax.axhline(y=-74.7360299558623, color='orange', label='PySCF mesh')
plt.legend()
plt.show()
