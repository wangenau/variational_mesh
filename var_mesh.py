import numpy as np
from pyscf import dft, gto
from pyscf.dft import radi


def mesh_error(coords, weights, atom_coord, atom_radius, alphas):
    '''Get the maximum mesh error for a few test functions.'''
    # shift atom position to the origin
    coords = coords - atom_coord
    mag = np.linalg.norm(coords, axis=1)
    # TODO: does not work properly somehow
    # mask = mag < atom_radius
    # mag = mag[np.array(mask)]
    # coords = coords[np.array(mask)]
    # weights = weights[np.array(mask)]
    err_max = 0
    for ia in alphas:
        for n in range(3):
            # numerical integration of (x*y*z)^n * exp(-a*r^2)
            val_num = 0
            for i in range(len(coords)):
                val_num += np.prod(coords[i]**n) * np.exp(-ia * mag[i]**2) * weights[i]
            # analytic solutions if we integrate from -inf to inf
            if n == 0:
                val_ana = np.sqrt(np.pi / ia)**3
            if n == 2:
                val_ana = np.sqrt(np.pi / (4 * ia**3))**3
            if n % 2 != 0:
                val_ana = 0
            if abs(val_num - val_ana) > err_max:
                err_max = abs(val_num - val_ana)
    return err_max


def var_mesh(mol, error=1e-3):
    '''Get a grid level whose error is smaller than a preset.'''
    # create grid class
    mesh = dft.Grids(mol)
    for il in range(10):
        print('Grid level: {}'.format(il))
        # build grid for a specific grid level
        mesh.level = il
        mesh.build()
        err_max = 0
        # sum up maximum errors for every atom
        for ia in range(mol.natm):
            coords = mesh.coords
            weights =  mesh.weights
            atom_coord = mol.atom_coord(ia)
            radius = radi.BRAGG_RADII[gto.charge(mol.atom_symbol(ia))]
            alphas = mol.bas_exp(ia)
            if len(alphas) > 2:
                alphas = [min(alphas), max(alphas)]
            err_max += mesh_error(coords, weights, atom_coord, radius, alphas)
        print('Max. error: {}'.format(err_max))
        if err_max < error:
            print('Error condition met.')
            return mesh
    print('Couldn\'t met error condition'.format(error))
    return mesh


if __name__ == "__main__":
    h2o = gto.Mole()
    h2o.atom = [
        ['C', ( 0    , 0, 0)],
        ['O', (-1.162, 0, 0)],
        ['O', ( 1.162, 0, 0)]]
    h2o.build()
    g = var_mesh(mol=h2o, error=1e-3)
    print(g.coords.shape)
