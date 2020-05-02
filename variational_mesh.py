#!/usr/bin/env python
import helpers
import numpy as np
from pyscf import dft, gto


def center_of_mass(coords, weights=None):
    '''Calculate the center of mass for a list of points and their weights.'''
    if not isinstance(weights, (list, np.ndarray)):
        weights = [1] * len(coords)
    com = np.full(len(coords[0]), 0)
    for i in range(len(coords)):
        com += coords[i] * weights[i]
    return com / sum(weights)


def main():
    mol = gto.M(atom='He 0, 0, 0', basis='sto-3g')
    mf = dft.RKS(mol).set(xc='lda,pw')
    mf.grids.level = 0
    mf.kernel()

    # get coordinates and weights and plot them
    atoms = mol.atom_coords()
    grids = mf.grids.coords
    weights = mf.grids.weights
    helpers.plot_mesh_3d(atoms=atoms, grids=grids, weights=weights)
    helpers.plot_mesh_2d(atoms=atoms, grids=grids, weights=weights, plane='xy')


if __name__ == "__main__":
    main()
