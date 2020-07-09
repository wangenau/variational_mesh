#!/usr/bin/env python
'''
Generate meshes with a predefined numerical error of the initial density.
'''

import numpy
import sys
from pyscf import dft, gto
from pyscf.dft import radi
from pyscf.dft.gen_grid import _default_ang, _default_rad
from pyscf.lib import logger


def gen_mesh(mesh, thres=1e-7):
    '''Minimize grid points for a given error for the initial density.

    Args:
        mesh :
            Grids object

    Kwargs:
        thres : scalar
            Maximum error for test functions in percent.

    Returns:
        Grid object
    '''
    # initialize logger
    verbose = mesh.verbose
    log = logger.Logger(sys.stdout, verbose)
    # create calculcator to calculate the initial density
    mf = dft.RKS(mesh.mol)
    mf.max_cycle = 0
    mf.grids = mesh
    mf.verbose = 0
    mf.mol.verbose = 0
    mf.grids.verbose = 0
    # enter coarse grid search
    types = atom_type(mesh.mol)
    log.note('Start coarse grid search.')
    for il in range(10):
        mf.grids = build_mesh(mf.grids, types, [il] * len(types))
        mf.kernel()
        error = mesh_error(mf)
        log.debug('[%d/9] Error: %.5E', il, error)
        if error < thres:
            log.note('Error condition met.\nLevel: %d', il)
            break
    if error >= thres:
        log.warn('Couldn\'t met error condition.\n'
                 'Use the largest grid level instead.')
        mf.mol.verbose = verbose
        return mf.grids
    # enter fine grid search
    level = il
    combs = get_combs(mesh.mol, level)
    counter = 0
    log.note('Start fine grid search.')
    for ic in combs:
        mf.grids = build_mesh(mf.grids, types, ic)
        mf.kernel()
        error = mesh_error(mf)
        log.debug('[%d/%d] Error: %.5E', counter, len(combs)-1, error)
        if error < thres:
            log.note('Error condition met.\nLevels per atom type:')
            for i in range(len(types)):
                log.note('\'%s\': %s', types[i], ic[i])
            break
        counter += 1
    if error >= thres:
        log.note('Couldn\'t enhance the mesh anymore.\n'
                 'Use level %d for every atom type instead.', il)
        mf.grids = build_mesh(mf.grids, types, [level] * len(types))
    # restore original verbose level
    mf.mol.verbose = verbose
    return mf.grids


def get_combs(mol, level):
    '''Generate possible grid level combinations for a fine grid search.'''
    types = atom_type(mol)
    amounts = atom_amount(mol)
    combs = numpy.empty((0, len(types)), int)
    # generate every possible combination
    for i in range(len(types)):
        tmp = [list(range(level, 10))] * len(types)
        tmp[i] = list(range(0, level))
        tmp = numpy.array(numpy.meshgrid(*tmp)).T.reshape(-1, len(types))
        combs = numpy.vstack((combs, tmp))
    # calculate grid points per combination
    grids = []
    for ic in combs:
        tmp = 0
        for i in range(len(types)):
            symb = types[i]
            chg = gto.charge(types[i])
            tmp += amounts[symb] * _default_rad(chg, ic[i]) * _default_ang(chg, ic[i])
        grids.append(tmp)
    # calculate the upper grid point boundary
    grids_max = 0
    for i in range(len(types)):
        symb = types[i]
        chg = gto.charge(symb)
        grids_max += amounts[symb] * _default_rad(chg, level) * _default_ang(chg, level)
    grids = numpy.asarray(grids)
    # remove combinations with more grids than the upper boundary
    mask = grids < grids_max
    grids = grids[mask]
    combs = combs[mask]
    # sort combinations with respect to grid points
    idx = numpy.argsort(grids)
    return combs[idx]


def build_mesh(mesh, types, levels):
    '''Build a mesh for given grid levels.'''
    for i in range(len(types)):
        symb = types[i]
        chg = gto.charge(symb)
        n_rad = _default_rad(chg, levels[i])
        n_ang = _default_ang(chg, levels[i])
        mesh.atom_grid[symb] = (n_rad, n_ang)
    return mesh.build()


def mesh_error(mf):
    '''Calculate the density error on a grid.'''
    mol = mf.mol
    dm = mf.make_rdm1()
    rho = mf._numint.get_rho(mol, dm, mf.grids, mf.max_memory)
    n = numpy.dot(rho, mf.grids.weights)
    return abs(mol.nelectron - n) / mol.nelectron


def atom_type(mol):
    '''Get types of atoms in the molecule.'''
    types = set()
    for ia in range(mol.natm):
        symb = mol.atom_symbol(ia)
        types.add(symb)
    return sorted(types)


def atom_amount(mol):
    '''Get amount of atoms in the molecule.'''
    amounts = {}
    for ia in range(mol.natm):
        symb = mol.atom_symbol(ia)
        amounts[symb] = sum(i.count(symb) for i in mol.atom)
    return amounts


if __name__ == "__main__":
    mol = gto.M(atom='C 0 0 0; O -1.162 0 0; O 1.162 0 0')
    mol.verbose = 9
    mesh = dft.Grids(mol)
    mesh = gen_mesh(mesh=mesh, thres=1e-7)
    print(mesh.coords.shape)
