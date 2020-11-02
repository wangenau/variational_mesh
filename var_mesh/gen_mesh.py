#!/usr/bin/env python
'''
Generate meshes by a numerical error threshold of the initial density.
'''

import numpy
import sys
from pyscf import dft, gto
from pyscf.lib import logger

# Angular grids that have to be used
ang_grids = dft.gen_grid.LEBEDEV_NGRID
rad = None
ang = None


def get_rad(symb, level):
    '''Get radial grids.'''
    try:
        return rad[symb][level]
    except:
        return dft.gen_grid._default_rad(gto.charge(symb), level)


def get_ang(symb, level):
    '''Get angular grids.'''
    try:
        return ang[symb][level]
    except:
        return dft.gen_grid._default_ang(gto.charge(symb), level)


def var_mesh(mesh, thres=1e-6, precise=True):
    '''Minimize grid points for a given error threshold for the initial density.

    Args:
        mesh :
            Grids object

    Kwargs:
        thres : scalar
            Maximum error of the initial density (normalized at one electron).

        precise : bool
            Use a precise search to find a better mesh that fits the threshold.
            Turned on by default (slower).

    Returns:
        Grids object
    '''
    # Initialize logger
    verbose = mesh.verbose
    log = logger.Logger(sys.stdout, verbose)
    # Create calculator to generate the initial density
    mf = dft.RKS(mesh.mol)
    mf.max_cycle = 0
    mf.grids = mesh
    # Silence other loggers
    mf.verbose = 0
    mf.mol.verbose = 0
    mf.grids.verbose = 0

    # Enter coarse grid search
    types = atom_types(mesh.mol)
    steps = get_steps()
    log.info('Start coarse grid search.')
    # Go through grid and raise atom grids for every atom
    for i in range(steps):
        mf.grids = build_mesh(mf.grids, types, [i] * len(types))
        mf.kernel()
        error = mesh_error(mf)
        log.debug('[%d/%d] Error = %.5e', i + 1, steps, error)
        if error < thres:
            log.info('Error condition met.')
            log.debug('Level = %d', i)
            break

    if error >= thres:
        log.warn('Couldn\'t met error condition.\n'
                 'Use the largest grid level instead.')
    elif precise:
        # Enter fine grid search
        level = i
        combs = get_combs(mesh.mol, level)
        steps = len(combs)
        log.info('Start fine grid search.')
        # Go through every possible grid level combination
        for i in range(steps):
            mf.grids = build_mesh(mf.grids, types, combs[i])
            mf.kernel()
            error = mesh_error(mf)
            log.debug('[%d/%d] Error = %.5e', i + 1, steps, error)
            if error < thres:
                log.info('Error condition met.')
                log.debug('Levels per atom type:')
                for j in range(len(types)):
                    log.debug('\'%s\' = %s', types[j], combs[i][j])
                break
        if error >= thres:
            log.info('Couldn\'t enhance the mesh anymore.')
            log.debug('Use level %d for every atom type instead.', level)
            mf.grids = build_mesh(mf.grids, types, [level] * len(types))

    log.note('Atom grids = %s', mf.grids.atom_grid)
    # Restore original verbose level
    mf.grids.verbose = verbose
    return mf.grids


def get_steps():
    '''Calculate the maximum number of optimization steps.'''
    if rad is None and ang is None:
        return 10
    elif rad is None:
        len_ang = min(len(i) for i in ang.values())
        return len_ang if len_ang <= 10 else 10
    elif ang is None:
        len_rad = min(len(i) for i in rad.values())
        return len_rad if len_rad <= 10 else 10
    else:
        len_rad = min(len(i) for i in rad.values())
        len_ang = min(len(i) for i in ang.values())
        return len_rad if len_rad <= len_ang else len_ang


def get_combs(mol, level):
    '''Generate possible grid level combinations for a fine grid search.'''
    steps = get_steps()
    types = atom_types(mol)
    amount = atom_amount(mol)
    combs = numpy.empty((0, len(types)), int)
    # Generate every possible combination
    for i in range(len(types)):
        tmp = [list(range(level, steps))] * len(types)
        tmp[i] = list(range(0, level))
        tmp = numpy.array(numpy.meshgrid(*tmp)).T.reshape(-1, len(types))
        combs = numpy.vstack((combs, tmp))
    # Calculate grid points per combination
    grids = []
    for ic in combs:
        tmp = 0
        for i in range(len(types)):
            symb = types[i]
            tmp += amount[symb] * get_rad(symb, ic[i]) * get_ang(symb, ic[i])
        grids.append(tmp)
    # Calculate the upper grid point boundary
    grids_max = 0
    for i in range(len(types)):
        symb = types[i]
        grids_max += amount[symb] * get_rad(symb, level) * get_ang(symb, level)
    grids = numpy.asarray(grids)
    # Remove combinations with more grids than the upper boundary
    mask = grids < grids_max
    grids = grids[mask]
    combs = combs[mask]
    # Sort combinations with respect to grid points
    idx = numpy.argsort(grids)
    return combs[idx]


def build_mesh(mesh, types, levels):
    '''Build a mesh for given grid levels.'''
    for i in range(len(types)):
        symb = types[i]
        n_rad = get_rad(symb, levels[i])
        n_ang = get_ang(symb, levels[i])
        mesh.atom_grid[symb] = (n_rad, n_ang)
    return mesh.build()


def mesh_error(mf):
    '''Calculate the density error per electron on a mesh.'''
    mol = mf.mol
    dm = mf.make_rdm1()
    # Account for different density matrix formats
    if dm.ndim == 3:
        dm = dm[0]
    rho = mf._numint.get_rho(mol, dm, mf.grids, mf.max_memory)
    n = numpy.dot(rho, mf.grids.weights)
    return abs(mol.nelectron - n) / mol.nelectron


def atom_types(mol):
    '''Get types of atoms in a molecule.'''
    types = set()
    for ia in range(mol.natm):
        symb = mol.atom_symbol(ia)
        types.add(symb)
    return sorted(types)


def atom_amount(mol):
    '''Get amount of atoms in a molecule.'''
    amount = {}
    for ia in range(mol.natm):
        symb = mol.atom_symbol(ia)
        amount[symb] = sum(i.count(symb) for i in mol.atom)
    return amount
