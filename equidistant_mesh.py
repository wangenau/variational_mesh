#!/usr/bin/env python
import helpers
import numpy as np
from pyscf import dft, gto
from pyscf.dft import radi


def gen_partition_equidist(mol, atom_grids_tab, radii_adjust=None,
                           atomic_radii=radi.BRAGG_RADII, becke_scheme=None):
    '''
    Generate an equidistant mesh around every atom of a given molecule.
    The mesh will have the form of a cube with edge lengths two times the
    atomic bragg radii.
    Signature is taken from pyscf.dft.gen_grid.gen_partition.
    '''
    nsteps = 10  # mesh points per axis
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


def main():
    mol = gto.M(atom='O 0, 0, 0; H 0, 1, 0; H 0, 0, 1', basis='sto-3g')
    mf = dft.RKS(mol).set(xc='lda,pw')
    mf.grids.gen_partition = gen_partition_equidist
    mf.kernel()

    # get coordinates of atoms and grid points and plot them
    atoms = mol.atom_coords()
    grids = mf.grids.coords
    helpers.plot_mesh_3d(atoms=atoms, grids=grids)
    helpers.plot_mesh_2d(atoms=atoms, grids=grids, plane='yz')


if __name__ == "__main__":
    main()
