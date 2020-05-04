#!/usr/bin/env python
import helpers
import numpy as np
from pyscf import dft, gto
from pyscf.lib.parameters import BOHR


def universal_box(mol, d=5):
    '''Calculate the universal box for a given molecule.

    Args:
        mol :
            :class:`Mole` object

    Kwargs:
        d : scalar
            Length that is added to both sides of the box in angstrom.

    Returns:
        Array with four coordinates that span the universal box.
    '''
    # calculate extra length for the box
    d = 2 * d / BOHR
    # calculate maximum length vectors and add d
    atm_coords = mol.atom_coords()
    xmax = [abs(max(atm_coords[:, 0]) - min(atm_coords[:, 0])) + d, 0, 0]
    ymax = [0, abs(max(atm_coords[:, 1]) - min(atm_coords[:, 1])) + d, 0]
    zmax = [0, 0, abs(max(atm_coords[:, 2]) - min(atm_coords[:, 2])) + d]
    # calculate geometric center of mass of the molecule
    atm_com = center_of_mass(atm_coords)
    # calculate center of mass of the box
    box_com = 0.5 * np.sum([xmax, ymax, zmax], 0)
    # create universal box and shift center of mass to the molecules one
    box = np.vstack(([0, 0, 0], xmax, ymax, zmax)) + (atm_com - box_com)
    return box


def center_of_mass(coords, weights=None):
    '''Calculate the center of mass for a list of points and their weights.'''
    if not isinstance(weights, (list, np.ndarray)):
        weights = [1] * len(coords)
    com = np.full(len(coords[0]), 0)
    for i in range(len(coords)):
        com = com + coords[i] * weights[i]
    return com / sum(weights)


def main():
    mol = gto.M(atom='O 0, 0, 0; H 0, 1, 0; H 0, 0, 1', basis='sto-3g')
    mf = dft.RKS(mol).set(xc='lda,pw')
    mf.grids.level = 0
    mf.kernel()

    atoms = mol.atom_coords()
    grids = mf.grids.coords
    cube = universal_box(mol)
    helpers.plot_mesh_3d(atoms=atoms, grids=grids, cubes=cube)
    helpers.plot_mesh_2d(atoms=atoms, grids=grids, cubes=cube, plane='yz')


if __name__ == "__main__":
    main()
