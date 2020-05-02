#!/usr/bin/env python
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
    mol = gto.M(atom='O 0, 0, 0; H 0, 1, 0; H 0, 0, 1', basis='sto-3g')
    atm_coords = np.asarray(mol.atom_coords())
    atm_masses = np.asarray(mol.atom_mass_list())
    atm_com = center_of_mass(atm_coords, atm_masses)
    print(atm_com)


if __name__ == "__main__":
    main()
