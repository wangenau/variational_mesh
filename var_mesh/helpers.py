#!/usr/bin/env python
'''
This file contains functions that may help to visualize, debug, or test
different meshes, but are not necessary to generate them.
'''

from mpl_toolkits.mplot3d import Axes3D
from pyscf.dft.radi import BRAGG_RADII
from pyscf.gto import charge
from var_mesh.gen_mesh import atom_type, get_combs
import matplotlib.pyplot as plt
import numpy as np

# Source: https://en.wikipedia.org/wiki/CPK_coloring
cpk_colors = {
    'H': 'w',
    'C': 'k',
    'N': 'b',
    'O': 'r',
    'S': 'y'
}


def plot_mesh_3d(mesh, weight=True, **kwargs):
    '''Plot atoms, grid points and parallelpipeds in a 3d plot.

    Args:
        mesh :
            Grids object

    Kwargs:
        weight : boolean
            Scale grid points by their weights if true.
    '''
    mol = mesh.mol
    coords = mesh.coords
    weights = 2
    if weight:
        weights = mesh.weights
    # generate information to plot the atoms
    atoms = mol.atom_coords()
    colors = []
    radii = []
    for ia in range(mol.natm):
        colors.append(cpk_colors.get(mol.atom_symbol(ia), 'magenta'))
        radii.append(BRAGG_RADII[charge(mol.atom_symbol(ia))] * 100)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=Axes3D.name)
    # plot atoms
    ax.scatter(atoms[:, 0], atoms[:, 1], atoms[:, 2], s=radii, c=colors,
               edgecolors='k', depthshade=0)
    # plot mesh points (by their weights) if available
    if isinstance(coords, np.ndarray):
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=weights, c='g')
    ax.set_xlabel('x axis')
    ax.set_ylabel('y axis')
    ax.set_zlabel('z axis')
    ax.autoscale()
    plt.show()
    return


def plot_mesh_2d(mesh, weight=True, plane='xy'):
    '''Project atoms, grid points and parallelpipeds to a given plane.

    Args:
        mesh :
            Grids object

    Kwargs:
        weight : boolean
            Scale grid points by their weights if true.

        plane : string
            Contains the plane to project on to.
    '''
    mol = mesh.mol
    coords = mesh.coords
    weights = 2
    if weight:
        weights = mesh.weights
    # generate information to plot the atoms
    atoms = mol.atom_coords()
    colors = []
    radii = []
    for ia in range(mol.natm):
        colors.append(cpk_colors.get(mol.atom_symbol(ia), 'magenta'))
        radii.append(BRAGG_RADII[charge(mol.atom_symbol(ia))] * 100)
    # dictionary to map input axes to their coordinates
    ax = {
        'x': 0,
        'y': 1,
        'z': 2
    }
    ax1, ax2 = list(plane.lower())
    # project atoms
    plt.scatter(atoms[:, ax[ax1]], atoms[:, ax[ax2]], s=radii, c=colors,
                edgecolors='k')
    # project mesh points (by their weights) if available
    if isinstance(coords, np.ndarray):
        plt.scatter(coords[:, ax[ax1]], coords[:, ax[ax2]], s=weights, c='g')
    plt.xlabel(ax1 + ' axis')
    plt.ylabel(ax2 + ' axis')
    plt.show()
    return


def plot_combs(mol, level):
    '''Display combinations for a fine grid search
       (assuming there are ten grid levels).

    Args:
        mol :
            Mole object

        level : scalar
            Contains the matching level of a coarse grid search.
    '''
    types = atom_type(mol)
    combs = get_combs(mol, level)
    if len(types) == 2:
        plt.scatter(level, level, c='r')
        plt.scatter(combs[:, 0], combs[:, 1])
        plt.xlim(-0.5, 9.5)
        plt.ylim(-0.5, 9.5)
        plt.xlabel('%s grid level' % types[0])
        plt.ylabel('%s grid level' % types[1])
    elif len(types) == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection=Axes3D.name)
        ax.scatter(level, level, level, c='r')
        ax.scatter(combs[:, 0], combs[:, 1], combs[:, 2])
        ax.set_xlim3d(-0.5, 9.5)
        ax.set_ylim3d(-0.5, 9.5)
        ax.set_zlim3d(-0.5, 9.5)
        ax.set_xlabel('%s grid level' % types[0])
        ax.set_ylabel('%s grid level' % types[1])
        ax.set_zlabel('%s grid level' % types[2])
    plt.show()
    return
