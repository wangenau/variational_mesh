#!/usr/bin/env python
'''
This module contains functions that are necessary to generate and optimize
meshes by numerical error thresholds of the initial-guess density.
'''

import numpy as np
from pyscf import dft, gto
from pyscf.lib import logger
from sys import stdout

# Lebedev angular grids that have to be used
ang_grids = dft.gen_grid.LEBEDEV_NGRID
rad = None
ang = None


def var_mesh(mesh, thres=1e-6, precise=True, mode='pyscf'):
    '''Minimize grid points for a given error threshold for the initial density.

    Args:
        mesh :
            Object of class :class:`Grids`.

    Kwargs:
        thres : float
            Maximum error of the initial-guess density (normalized to one
            electron).

        precise : bool
            Disable the coarse grid search when set to ``False``. This will
            speed up the mesh generation but results in more grid points.

        mode : str
            Set the output format. If the code does not support different grids
            for different atom types, only the largest atom grid will be used
            in an optimization step. Can be one of ``'pyscf'``, ``'erkale'`` or
            ``'gamess'``.

    Returns:
        Object of class :class:`Grids`.
    '''
    # Initialize logger
    verbose = mesh.verbose
    log = logger.Logger(stdout, verbose)
    # Handle user inputs
    if not isinstance(mesh, dft.gen_grid.Grids):
        log.error('The mesh parameter has to be a Grids object.')
    if not isinstance(thres, float) and not isinstance(thres, int):
        log.error('The threshold parameter has to be a number.')
    if thres < 0:
        log.error('The threshold parameter has to be positive.')
    if not isinstance(precise, bool):
        log.error('The precise parameter has to be a boolean.')
    if not isinstance(mode, str):
        log.error('The mode parameter has to be a string.')
    mode = mode.lower()
    supported = ('pyscf', 'erkale', 'gamess')
    if mode not in supported:
        log.error('The mode parameter has to be one of \'%s\'.',
                  '\', \''.join(supported))
    if mode != 'pyscf' and precise:
        log.warn('The precise parameter has no effect when using \'%s\' mode.',
                 mode)

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
        mf.grids = build_mesh(mf.grids, types, [i] * len(types), mode)
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
    elif precise and mode == 'pyscf':
        # Enter fine grid search
        level = i
        combs = get_combs(mesh.mol, level)
        steps = len(combs)
        log.info('Start fine grid search.')
        # Go through every possible grid level combination
        for i in range(steps):
            mf.grids = build_mesh(mf.grids, types, combs[i], mode)
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
            mf.grids = build_mesh(mf.grids, types, [level] * len(types), mode)
    if mode == 'erkale':
        grid = list(mf.grids.atom_grid.values())[0]
        log.note('ERKALE grid: DFTGrid %d -%d', grid[0], grid[1])
    elif mode == 'gamess':
        grid = list(mf.grids.atom_grid.values())[0]
        log.note('GAMESS grid: $dft nrad=%d nleb=%d $end', grid[0], grid[1])
    else:
        log.note('PySCF grid: atom_grid = %s', mf.grids.atom_grid)
    # Restore original verbose level
    mf.grids.verbose = verbose
    return mf.grids


def get_rad(symb, level):
    '''Get radial grids.

    Args:
        symb : str
            Atom type identifier.

        level : int
            Grid level of atom type.

    Returns: int
        Amount of radial grids.
    '''
    try:
        return rad[symb][level]
    except (KeyError, TypeError):
        return dft.gen_grid._default_rad(gto.charge(symb), level)


def get_ang(symb, level):
    '''Get angular grids.

    Args:
        symb : str
            Atom type identifier.

        level : int
            Grid level of atom type.

    Returns: int
        Amount of angular grids.
    '''
    try:
        return ang[symb][level]
    except (KeyError, TypeError):
        return dft.gen_grid._default_ang(gto.charge(symb), level)


def get_steps():
    '''Get the maximum number of optimization steps.

    Returns: int
        Largest grid level.
    '''
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
    '''Generate possible grid level combinations for a fine grid search.

    Args:
        mol :
            Object of class :class:`Mole`.

        level : int
            Grid level of the found match in the coarse grid search.

    Returns: ndarray
        Possible combinations of grid levels.
    '''
    steps = get_steps()
    types = atom_types(mol)
    amount = atom_amount(mol)
    combs = np.empty((0, len(types)), int)
    # Generate every possible combination
    for i in range(len(types)):
        tmp = [list(range(level, steps))] * len(types)
        tmp[i] = list(range(0, level))
        tmp = np.array(np.meshgrid(*tmp)).T.reshape(-1, len(types))
        combs = np.vstack((combs, tmp))
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
    grids = np.asarray(grids)
    # Remove combinations with more grids than the upper boundary
    mask = grids < grids_max
    grids = grids[mask]
    combs = combs[mask]
    # Sort combinations with respect to grid points
    idx = np.argsort(grids)
    return combs[idx]


def build_mesh(mesh, types, levels, mode):
    '''Build mesh for given grid levels.

    Args:
        mesh :
            Object of class :class:`Grids`.

        types : list
             Atom type identifiers.

        levels : list
            Grid levels per atom type.

        mode : str
            Set the output format. If the code does not support different grids
            for different atom types, only the largest atom grid will be used
            in an optimization step. Can be one of ``'pyscf'``, ``'erkale'`` or
            ``'gamess'``.

    Returns:
        Object of class :class:`Grids`.
    '''
    if mode == 'pyscf':
        for i in range(len(types)):
            symb = types[i]
            n_rad = get_rad(symb, levels[i])
            n_ang = get_ang(symb, levels[i])
            mesh.atom_grid[symb] = (n_rad, n_ang)
    else:
        n_rads = []
        n_angs = []
        for i in range(len(types)):
            symb = types[i]
            n_rads.append(get_rad(symb, levels[i]))
            n_angs.append(get_ang(symb, levels[i]))
        # Get the index for the largest grid and use it for every atom type
        idx = np.argmax(np.multiply(n_rads, n_angs))
        for it in types:
            mesh.atom_grid[it] = (n_rads[idx], n_angs[idx])
    return mesh.build()


def mesh_error(mf):
    '''Calculate the error per electron when integrating the electron density.

    Args:
        mol :
            Object of class :class:`RKS` or :class:`UKS`.

    Returns: float
        Mesh error.
    '''
    mol = mf.mol
    dm = mf.make_rdm1()
    # Density matrix for open-shell systems
    if dm.ndim != 2:
        dm = dm[0] + dm[1]
    rho = mf._numint.get_rho(mol, dm, mf.grids, mf.max_memory)
    n = np.dot(rho, mf.grids.weights)
    return abs(mol.nelectron - n) / mol.nelectron


def atom_types(mol):
    '''Get types of atoms in a molecule.

    Args:
        mol :
            Object of class :class:`Mole`.

    Returns: set
        Unique atom type identifiers.
    '''
    types = set()
    for ia in range(mol.natm):
        symb = mol.atom_symbol(ia)
        types.add(symb)
    return sorted(types)


def atom_amount(mol):
    '''Get the amount of atoms in a molecule.

    Args:
        mol :
            Object of class :class:`Mole`.

    Returns: dict
        Atom type identifiers as keys and amounts as values.
    '''
    amount = {}
    for ia in range(mol.natm):
        symb = mol.atom_symbol(ia)
        amount[symb] = sum(i[0].count(symb) for i in mol.atom)
    return amount
