import itertools
import numpy as np
from pyscf import dft, gto
from pyscf.dft.gen_grid import gen_atomic_grids, _default_ang, _default_rad


def mesh_error(coords, weights, atom_coord, alphas):
    '''Get the maximum mesh error for a few test functions.'''
    # analytic solutions if we integrate (x*y*z)^n * exp(-a*r^2) dxdydz
    # from -inf to inf for a specific n
    solutions = {
        0: lambda ia: np.sqrt(np.pi / ia)**3,
        1: lambda ia: 0,
        2: lambda ia: np.sqrt(np.pi / (4 * ia**3))**3
    }

    # shift atom position to the origin
    coords = coords - atom_coord
    err_max = 0
    for ia in alphas:
        for n in range(3):
            # numerical integration of (x*y*z)^n * exp(-a*r^2)
            val_num = 0
            for i in range(len(coords)):
                val_num += np.prod(coords[i]**n) * \
                    np.exp(-ia * np.linalg.norm(coords[i])**2) * weights[i]
            val_ana = solutions[n](ia)
            # select the highest error
            if abs(val_num - val_ana) > err_max:
                err_max = abs(val_num - val_ana)
    return err_max


def var_mesh(mol, error=1e-3):
    '''Get a grid level whose error is smaller than a preset.'''
    # list of atom types
    atoms = []
    for ia in range(len(mol.atom)):
        key = mol.atom_symbol(ia)
        atoms.append(key)
    atoms = list(set(atoms))
    print('List of atoms: {}\n'.format(atoms))

    # save the amount of atoms in a dict
    atom_amount = {}
    for key in atoms:
        amount = sum(ia.count(key) for ia in mol.atom)
        atom_amount[key] = amount
    print('Amount of atoms: {}\n'.format(atom_amount))

    # create an initial mesh grid
    mesh = dft.Grids(mol)
    mesh.level = 0
    mesh.build()
    coords = mesh.coords
    weights = mesh.weights

    # save the index for every atom type with the maximum error in a dict
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
    print('Index for atom type: {}\n'.format(atom_index))

    # lists that contain the error/#grid points per atom type and grid level
    errors = []
    grids = []
    # start loop to go through every grid level
    for il in range(10):
        print('Grid level: {}'.format(il))
        err_max = 0
        tmp_errors = []
        tmp_grids = []
        # build grid for a specific grid level
        mesh.level = il
        mesh.build()
        coords = mesh.coords
        weights = mesh.weights
        # generate atomic grids to get the number of grid points
        mesh_dict = gen_atomic_grids(mol, level=il)
        for key in atoms:
            atom_coord = mol.atom_coord(atom_index[key])
            alphas = mol.bas_exp(atom_index[key])
            if len(alphas) > 2:
                alphas = [min(alphas), max(alphas)]
            err = mesh_error(coords, weights, atom_coord, alphas)
            # add error to max error, multiplied with the amount of atoms
            err *= atom_amount[key]
            err_max += err
            # append errors and #grid points per atom type
            tmp_errors.append(err)
            tmp_grids.append(atom_amount[key] * len(mesh_dict[key][0]))
        errors.append(tmp_errors)
        grids.append(tmp_grids)
        print('Max. error: {}'.format(err_max))
        if err_max < error:
            print('Error condition met.\n')
            break
    if il == 9 and err_max > error:
        print('Couldn\'t met error condition.\n')
        return mesh

    errors = np.asarray(errors)
    grids = np.asarray(grids)
    # generate every possibility to combine grids
    combinations = list(itertools.product(*[range(il + 1)] * len(atoms)))
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
    for ia in range(len(atoms)):
        n_rad = _default_rad(gto.charge(atoms[ia]), min_levels[ia])
        n_ang = _default_ang(gto.charge(atoms[ia]), min_levels[ia])
        mesh.atom_grid[atoms[ia]] = (n_rad, n_ang)
        print('Grid level for \'{}\': {}'.format(atoms[ia], min_levels[ia]))
    mesh.build()

    # final error check
    # coords = mesh.coords
    # weights = mesh.weights
    # err = 0
    # for ia in range(mol.natm):
    #     atom_coord = mol.atom_coord(ia)
    #     alphas = mol.bas_exp(ia)
    #     if len(alphas) > 2:
    #         alphas = [min(alphas), max(alphas)]
    #     err_max += mesh_error(coords, weights, atom_coord, alphas)
    # print('Max. final error: {}'.format(err_max))
    return mesh


if __name__ == "__main__":
    mol = gto.Mole()
    mol.atom = [
        ['C', ( 0    , 0, 0)],
        ['O', (-1.162, 0, 0)],
        ['O', ( 1.162, 0, 0)]]
    mol.build()
    g = var_mesh(mol=mol, error=1e-3)
    print(g.coords.shape)
