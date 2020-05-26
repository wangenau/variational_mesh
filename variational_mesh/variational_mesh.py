#!/usr/bin/env python
import numpy as np
import sys
from pyscf import dft, gto
from pyscf.dft import radi
from pyscf.dft.gen_grid import *
from pyscf.dft.gen_grid import _default_ang, _default_rad
from pyscf.lib import logger


def mesh_error(coords, weights, atom_coord, alphas):
    '''Calculate the maximum grid error for test functions.

    Args:
        coords : array
            Contains the coordinates of every grid point.

        weights : array
            Contains the weights of every grid point.

        atom_coord : array
            Contains the coordinates of the atom position.

        alphas : list
            Contains the coefficients of the atoms basis set.

    Returns:
        Maximum error in percent from test function integrations.
    '''
    # analytic solutions for the integration of exp(-a*r^2)*(x*y*z)^n from
    # -inf to inf for different n
    solutions = {
        0: lambda ia: np.sqrt(np.pi / ia)**3,
        1: lambda ia: 0,
        2: lambda ia: np.sqrt(np.pi / (4 * ia**3))**3
    }

    # shift atom position to the origin
    coords = coords - atom_coord
    err_max = 0
    for ia in alphas:
        for n in range(len(solutions)):
            # numerical integrate exp(-a*r^2)*(x*y*z)^n
            val_num = 0
            for i in range(len(coords)):
                val_num += np.prod(coords[i]**n) * \
                    np.exp(-ia * np.linalg.norm(coords[i])**2) * weights[i]
            # select the highest error (in percent)
            val_ana = solutions[n](ia)
            err = abs(val_num - val_ana)
            if err > err_max:
                err_max = err
    return err_max


def variational_mesh(mol, error=1e-3, radi_method=radi.treutler,
                     prune=nwchem_prune):
    '''Minimize grid point amount for a given maximum error for test functions.

    Args:
        mol :
            Mole object

    Kwargs:
        error : scalar
            Maximum error for test functions in percent.

        radi_method : function(n) => (rad_grids, rad_weights)
            scheme for radial grids, can be one of
            | radi.treutler  (default)
            | radi.delley
            | radi.mura_knowles
            | radi.gauss_chebyshev

        prune : function(nuc, rad_grids, n_ang) => list_n_ang_for_each_rad_grid
            scheme to reduce number of grids, can be one of
            | gen_grid.nwchem_prune  (default)
            | gen_grid.sg1_prune
            | gen_grid.treutler_prune
            | None : to switch off grid pruning

    Returns:
        Grid coordinates and weights arrays.
    '''
    # initialize logger and shut up other loggers
    log = logger.Logger(sys.stdout, mol.verbose)
    mol.verbose = 0

    # list of individual atoms
    atoms = []
    for ia in range(len(mol.atom)):
        key = mol.atom_symbol(ia)
        atoms.append(key)
    atoms = list(set(atoms))
    log.info('List of atoms: %s', atoms)

    # save the amount of atoms in a dict
    atom_amount = {}
    for key in atoms:
        amount = sum(ia.count(key) for ia in mol.atom)
        atom_amount[key] = amount
    log.info('Amount of atoms: %s', atom_amount)

    # create an initial grid at the lowest grid level
    mesh = dft.Grids(mol)
    mesh.level = 0
    mesh.radi_method = radi_method
    mesh.prune = prune
    mesh.build()
    coords = mesh.coords
    weights = mesh.weights
    mesh_dict = gen_atomic_grids(mol, level=0, radi_method=radi_method, prune=prune)

    # lists that contain the error/grid point number per atom type and grid level
    errors = [[]]
    grids = [[]]
    # save the index for every atom type with the maximum error in a dict
    # also run the first grid level explicitly and save the results
    atom_index = {}
    for key in atoms:
        err_max = 0
        for ia in range(mol.natm):
            if key == mol.atom_symbol(ia):
                atom_coord = mol.atom_coord(ia)
                alphas = mol.bas_exp(ia)
                if len(alphas) > 2:
                    alphas = [min(alphas), max(alphas)]
                err = mesh_error(coords, weights, atom_coord, alphas)
                if err > err_max:
                    err_max = err
                    atom_index[key] = ia
        errors[0].append(atom_amount[key] * err_max)
        grids[0].append(atom_amount[key] * len(mesh_dict[key][0]))
    log.debug('Index per atom type: %s', atom_index)
    log.debug('Grid level: 0')
    log.debug('Max. error: %f', sum(errors[0]))

    # go through every grid level until the error condition is met
    for il in range(1, 10):
        log.debug('Grid level: %d', il)
        err_max = 0
        tmp_errors = []
        tmp_grids = []
        # build grid for a specific grid level
        mesh.level = il
        mesh.build()
        coords = mesh.coords
        weights = mesh.weights
        # generate atomic grids to get the number of grid points
        mesh_dict = gen_atomic_grids(mol, level=il, radi_method=radi_method, prune=prune)
        for key in atoms:
            atom_coord = mol.atom_coord(atom_index[key])
            alphas = mol.bas_exp(atom_index[key])
            if len(alphas) > 2:
                alphas = [min(alphas), max(alphas)]
            err = mesh_error(coords, weights, atom_coord, alphas)
            # add error to max error, multiplied with the amount of atoms
            err *= atom_amount[key]
            err_max += err
            # append errors and  grid point amount per atom type
            tmp_errors.append(err)
            tmp_grids.append(atom_amount[key] * len(mesh_dict[key][0]))
        errors.append(tmp_errors)
        grids.append(tmp_grids)
        log.debug('Max. error: %f', err_max)
        if err_max < error:
            log.info('Error condition met.')
            break
    if il == 9 and err_max > error:
        log.warn('Couldn\'t met error condition.')
        return mesh

    errors = np.asarray(errors)
    log.debug('Errors per atom and grid level:\n{}'.format(errors))
    grids = np.asarray(grids)
    log.debug('Grid point amount per atom and grid level:\n{}'.format(grids))
    # generate every possible combimiation to create grids
    combinations = np.array(np.meshgrid(*[range(il + 1)] * len(atoms))).T.reshape(-1, len(atoms))
    min_levels = [il] * len(atoms)
    min_grids = sum(grids[-1])
    # search for the smallest amount of grid points that is under our error
    for ic in combinations:
        tmp_grids = sum(grids[ic, range(len(atoms))])
        tmp_error = sum(errors[ic, range(len(atoms))])
        if tmp_grids < min_grids and tmp_error < error:
            min_levels = ic
            min_grids = tmp_grids

    # build the final grid with the found atomic grid levels
    mesh = dft.Grids(mol)
    mesh.radi_method = radi_method
    mesh.prune = prune
    for ia in range(len(atoms)):
        n_rad = _default_rad(gto.charge(atoms[ia]), min_levels[ia])
        n_ang = _default_ang(gto.charge(atoms[ia]), min_levels[ia])
        mesh.atom_grid[atoms[ia]] = (n_rad, n_ang)
        log.info('Grid level for \'%s\': %d', atoms[ia], min_levels[ia])
    mesh.build()
    return mesh


if __name__ == "__main__":
    mol = gto.Mole()
    mol.atom = [
        ['C', ( 0    , 0, 0)],
        ['O', (-1.162, 0, 0)],
        ['O', ( 1.162, 0, 0)]]
    mol.build()
    mol.verbose = 5
    mesh = variational_mesh(mol=mol, error=1e-3)
    print(mesh.coords.shape)
