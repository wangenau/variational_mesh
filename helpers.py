#!/usr/bin/env python
'''
This file contains functions that may help to visualize, debug, or test
different meshes, but are not necessary needed to generate them.
'''

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np


def plot_mesh_3d(atoms, grids=None, weights=None, cubes=None, **kwargs):
    '''Plot atoms, grid points and parallelpipeds in a 3d plot.

    Args:
        atoms : array
            Contains coordinates of atom positions.

    Kwargs:
        grids : array
            Contains coordinates of grid points.

        weights : array or list
            Contains weights of every grid point.

        cubes : array
            Contains coordinates of four points that span parallelpipeds.
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=Axes3D.name)
    # plot atoms
    ax.scatter(atoms[:, 0], atoms[:, 1], atoms[:, 2], s=50, c='r')
    # plot mesh points (by their weights) if available
    if isinstance(grids, np.ndarray):
        if not isinstance(weights, (list, np.ndarray)):
            weights = 1
        ax.scatter(grids[:, 0], grids[:, 1], grids[:, 2], s=weights, c='g')
    # plot cubes if available
    if isinstance(cubes, np.ndarray):
        # if only one cube available, add an extra dimension
        if cubes.ndim == 2:
            cubes = np.array([cubes])
        for ic in cubes:
            faces = hexagon_faces(ic)
            # do not call with color argument, it will override alpha channel
            cube = Poly3DCollection(faces, alpha=0.1)
            cube.set_color('b')
            ax.add_collection3d(cube)
    ax.set_xlabel('x axis')
    ax.set_ylabel('y axis')
    ax.set_zlabel('z axis')
    ax.autoscale()
    plt.show()
    return


def plot_mesh_2d(atoms, grids=None, weights=None, cubes=None, plane='xy'):
    '''Project atoms, grid points and parallelpipeds to a given plane.

    Args:
        atoms : array
            Contains coordinates of atom positions.

    Kwargs:
        grids : array
            Contains coordinates of grid points.

        weights : array or list
            Contains weights of every grid point.

        cubes : array
            Contains coordinates of four points that span parallelpipeds.

        plane : string
            Contains the plane to project on to.
    '''
    # dictionary to map input axes to their coordinates
    ax = {
        'x': 0,
        'y': 1,
        'z': 2
    }
    ax1, ax2 = list(plane.lower())
    # project atoms
    plt.scatter(atoms[:, ax[ax1]], atoms[:, ax[ax2]], s=50, c='r')
    # project mesh points (by their weights) if available
    if isinstance(grids, np.ndarray):
        if not isinstance(weights, (list, np.ndarray)):
            weights = 1
        plt.scatter(grids[:, ax[ax1]], grids[:, ax[ax2]], s=weights, c='g')
    # project cubes if available
    if isinstance(cubes, np.ndarray):
        # if only one cube available, add an extra dimension
        if cubes.ndim == 2:
            cubes = np.array([cubes])
        for ic in cubes:
            areas = hexagon_faces(ic)
            for ia in areas:
                area = np.vstack(ia)
                plt.fill(area[:, ax[ax1]], area[:, ax[ax2]], alpha=0.1, c='b')
    plt.xlabel(ax1 + ' axis')
    plt.ylabel(ax2 + ' axis')
    plt.show()
    return


def hexagon_faces(cube):
    '''Calculate corners of every face of an arbitrary parallelpiped.

    Args:
        cube : array or list
            Contains coordinates of four points that span a parallelpiped.

    Returns:
        Array with coordinates of four corner points for every face.
    '''
    points = np.asarray(cube)
    points = np.vstack((points,
                        points[1] + points[2] - points[0],
                        points[1] + points[3] - points[0],
                        points[2] + points[3] - points[0],
                        points[1] + points[2] + points[3] - 2 * points[0]))
    faces = [[points[0], points[2], points[4], points[1]],
             [points[0], points[3], points[5], points[1]],
             [points[1], points[5], points[7], points[4]],
             [points[2], points[6], points[3], points[0]],
             [points[3], points[6], points[7], points[5]],
             [points[4], points[2], points[6], points[7]]]
    return faces
