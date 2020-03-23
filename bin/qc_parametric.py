#!/usr/bin/env python
'''
Given a quadratic surface, and a FEM mesh, generate a .msh file containing
a visualization of the parameterized surface over the head model for QC

Usage:
    parameterize_surface_patch <qc_surf> <C> <R> <bounds> <out>

Arguments:
    <qc_surf>                             Coordinates for a surface patch
    <C>                                   Quadratic constants
    <R>                                   Affine matrix
    <bounds>                              Sampling bounds
    <out>                                 Output filename
'''

from docopt import docopt
from fieldopt import geolib as gl
import numpy as np
import gmsh


def main():

    args = docopt(__doc__)

    qc_surf = args['<qc_surf>']
    C = np.load(args['<C>'])
    R = np.load(args['<R>'])
    b = np.load(args['<bounds>'])
    out = args['<out>']

    # Construct parameteric surface coordinates
    X, Y = np.meshgrid(np.linspace(b[0, 0] + 1, b[0, 1] - 1),
                       np.linspace(b[1, 0] + 1, b[1, 1] - 1))
    XX = X.flatten()
    YY = Y.flatten()
    poly_arr = np.c_[np.ones(XX.shape[0]), XX, YY, XX * YY, XX * XX, YY * YY]
    Z = np.dot(poly_arr, C)
    poly_coords = np.matmul(R, np.c_[XX, YY, Z].T).T

    # Load in mesh
    gmsh.initialize()
    gmsh.open(qc_surf)
    entity = gmsh.model.getEntities()[0]
    n_tag, n_coord, params = gmsh.model.mesh.getNodes(entity[0], entity[1])

    samp_surf_coords = poly_coords.flatten(order='C')
    samp_surf_nodes = np.arange(max(n_tag), max(n_tag) + poly_coords.shape[0])
    grid = samp_surf_nodes.reshape((50, 50))
    n = grid.shape[0]

    # Construct mesh
    trig_list = np.zeros((6 * (n - 1)**2), dtype=np.int64)
    for j in np.arange(0, grid.shape[0] - 1):
        for i in np.arange(0, grid.shape[1] - 1):
            # Upper triangular
            trig_list[6 * i + 6 * (n - 1) * j] = grid[(j, i)]
            trig_list[6 * i + 1 + 6 * (n - 1) * j] = grid[(j, i + 1)]
            trig_list[6 * i + 2 + 6 * (n - 1) * j] = grid[(j + 1, i + 1)]

            # Lower triangular
            trig_list[6 * i + 3 + 6 * (n - 1) * j] = grid[(j + 1, i)]
            trig_list[6 * i + 4 + 6 * (n - 1) * j] = grid[(j, i)]
            trig_list[6 * i + 5 + 6 * (n - 1) * j] = grid[(j + 1, i + 1)]

    gmsh.initialize()
    gmsh.model.add('sampling_surf')
    tag = gmsh.model.addDiscreteEntity(2, 3002)
    gmsh.model.mesh.setNodes(2,
                             tag,
                             nodeTags=samp_surf_nodes,
                             coord=samp_surf_coords)
    gmsh.model.mesh.setElements(
        2,
        tag, [2],
        elementTags=[range(1,
                           len(trig_list) // 3 + 1)],
        nodeTags=[trig_list])
    gmsh.write(out)
    gmsh.finalize()


if __name__ == '__main__':
    main()
