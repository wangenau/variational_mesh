#!/usr/bin/env python
'''
This module contains functions to visualize different meshes, but are not
necessary to generate them.
'''

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pyscf.dft.gen_grid import Grids
from pyscf.dft.radi import BRAGG_RADII
from pyscf.gto import charge
from pyscf.lib import logger
from sys import stdout

# Color some atoms with their respective CPK color
# https://en.wikipedia.org/wiki/CPK_coloring
cpk_colors = {
    'H': 'w',
    'C': 'k',
    'N': 'b',
    'O': 'r',
    'S': 'y'
}


def plot_mesh_3d(mesh, weight=False, **kwargs):
    '''Plot atoms and grid points in a 3d plot.

    Args:
        mesh :
            Object of class :class:`Grids`.

    Kwargs:
        weight : bool or int
            Scale grid points by their weights when set to ``True``. Integers
            will scale grid points.
    '''
    # Initialize logger
    verbose = mesh.verbose
    log = logger.Logger(stdout, verbose)
    # Handle user inputs
    if not isinstance(mesh, Grids):
        log.error('The mesh parameter has to be a Grids object.')
    if not isinstance(weight, bool) and not isinstance(weight, int):
        log.error('The weight parameter has to be a boolean or an integer.')

    mol = mesh.mol
    coords = mesh.coords
    atoms = mol.atom_coords()
    if weight and isinstance(weight, bool):
        weights = abs(mesh.weights)
    elif weight and isinstance(weight, int):
        weights = weight
    else:
        weights = 2
    # Get the atom colors and radii (scale with 100 to make it look good)
    colors = []
    radii = []
    for ia in range(mol.natm):
        colors.append(cpk_colors.get(mol.atom_symbol(ia), 'magenta'))
        radii.append(BRAGG_RADII[charge(mol.atom_symbol(ia))] * 100)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=Axes3D.name)
    # Plot atoms
    ax.scatter(atoms[:, 0], atoms[:, 1], atoms[:, 2], s=radii, c=colors,
               edgecolors='k', depthshade=0)
    # Plot mesh points
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=weights, c='g')
    ax.set_xlabel('x-axis', fontsize=12)
    ax.set_ylabel('y-axis', fontsize=12)
    ax.set_zlabel('z-axis', fontsize=12)
    plt.tight_layout()
    plt.show()
    return


def plot_mesh_2d(mesh, weight=False, plane='xy'):
    '''Project atoms and grid points to a given plane in a 2d plot.

    Args:
        mesh :
            Object of class :class:`Grids`.

    Kwargs:
        weight : bool or int
            Scale grid points by their weights when set to ``True``. Integers
            will scale grid points.

        plane : str
            Contains the plane to project on to. Needs two elements of ``'x'``,
            ``'y'``, and ``'z'``.
    '''
    # Dictionary to map input axes to their coordinates
    ax = {
        'x': 0,
        'y': 1,
        'z': 2
    }

    # Initialize logger
    verbose = mesh.verbose
    log = logger.Logger(stdout, verbose)
    # Handle user inputs
    if not isinstance(mesh, Grids):
        log.error('The mesh parameter has to be a Grids object.')
    if not isinstance(weight, bool) and not isinstance(weight, int):
        log.error('The weight parameter has to be a boolean or an integer.')
    if not isinstance(plane, str):
        log.error('The plane parameter has to be a string.')
    if len(plane) != 2:
        log.error('The plane parameter has to be a string of length 2.')
    if not set(plane).issubset(ax.keys()):
        log.error('The plane parameter only accepts combinations of \'%s\'.',
                  '\', \''.join(ax.keys()))

    ax1, ax2 = list(plane.lower())
    mol = mesh.mol
    coords = mesh.coords
    atoms = mol.atom_coords()
    if weight and isinstance(weight, bool):
        weights = abs(mesh.weights)
    elif weight and isinstance(weight, int):
        weights = weight
    else:
        weights = 2
    # Get the atom colors and radii (scale with 100 to make it look good)
    colors = []
    radii = []
    for ia in range(mol.natm):
        colors.append(cpk_colors.get(mol.atom_symbol(ia), 'magenta'))
        radii.append(BRAGG_RADII[charge(mol.atom_symbol(ia))] * 100)
    # Project atoms
    plt.scatter(atoms[:, ax[ax1]], atoms[:, ax[ax2]], s=radii, c=colors,
                edgecolors='k')
    # Project mesh points
    plt.scatter(coords[:, ax[ax1]], coords[:, ax[ax2]], s=weights, c='g')
    plt.xlabel(ax1 + '-axis', fontsize=12)
    plt.ylabel(ax2 + '-axis', fontsize=12)
    plt.tight_layout()
    plt.show()
    return
