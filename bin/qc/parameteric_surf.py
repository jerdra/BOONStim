#!/usr/bin/env python
import argparse
import logging
from collections import namedtuple

import numpy as np
from fieldopt.geometry import mesh as gm
from fieldopt.geometry import geometry as gl
from fieldopt.geometry.mesh_wrapper import HeadModel
from fieldopt.geometry.domains import QuadraticDomain

import meshplot as mp

HEAD_ENTITIES = [(2, 5), (2, 1005)]
GM_ENTITIES = [(2, 2), (2, 1002)]

SurfaceMesh = namedtuple("SurfaceMesh", ["nodes", "coords", "triangles"])
DistResult = namedtuple("DistResult", ["source", "target", "distance"])

logging.basicConfig(format="%(asctime)s [BOONSTIM CORTEX2SCALP]:  %(message)s",
                    datefmt="%Y-%m-%d %I:%M:%S %p",
                    level=logging.INFO)

mp.website()


def decompose_gmsh(f_mesh, entity):
    '''
    Given a path to a mesh file load the:
        nodes
        coordinates
        triangles

    For the specified surface.

    Arguments:
        f_mesh          Path to GMSH mesh file
        entity          Tuple describing gmsh (dim, tag)
    '''

    # Load vertex data
    nodes, coords, _ = load_surf_verts(f_mesh, entity)
    trigs = load_surf_trigs(f_mesh, entity)

    # Re-normalize node/trig identities
    trigs = trigs - nodes.min()
    nodes = nodes - nodes.min()
    return SurfaceMesh(nodes, coords, trigs)


def load_surf_verts(f_msh, entities):
    '''
    surf:       Path to gmsh MSH file
    entities:   List of entities to attempt of form (dim, tag)
    '''

    for dim, tag in entities:

        try:
            nodes, coords, params = gm.load_gmsh_nodes(f_msh, (dim, tag))
        except ValueError:
            continue
        else:
            return nodes, coords, params

    logging.error("Could not properly load Mesh! Check entity tags!")
    raise ValueError


def load_surf_trigs(f_msh, entities):
    '''
    surf:       Path to gmsh MSH file
    entities:   List of entities to attempt of form (dim, tag)
    '''

    for dim, tag in entities:

        try:
            _, _, trigs = gm.load_gmsh_elems(f_msh, (dim, tag))
        except ValueError:
            continue
        else:
            return trigs.reshape(-1, 3)

    logging.error("Could not properly load Mesh! Check entity tags!")
    raise ValueError


def gen_mshplot_html(brain, texture, param_surf, out_html):
    '''
    Generate a visualization of the centroid to head projection problem

    Arguments:
       brain:           SurfaceMesh of brain
       param_surf:      Parameteric surface to visualize
       out_html:        Output HTML file
    '''

    p = mp.plot(brain.coords, brain.triangles, c=texture)

    p.add_mesh(param_surf.coords,
               param_surf.triangles,
               c=np.array([0.7, 0.7, 0.7]))
    p._Viewer__objects[1]['material'].transparent = True
    p._Viewer__objects[1]['material'].opacity = 0.5
    p._Viewer__objects[1]['material'].metalness = 0.1

    p.save(out_html)


def construct_parameteric_mesh(C, A, bounds, grid_res=50):
    '''
    Construct a Mesh representation of the parameteric surface
    used for continuous sampling

    Arguments:
        C           Quadratic surface constants
        A           Affine transformation into mesh space
        bounds      Limit bounds describing the extent of the surface

    Returns:
        SurfaceMesh A SurfaceMesh object for visualizing the parameteric mesh
    '''

    X, Y = np.meshgrid(np.linspace(bounds[0, 0], bounds[0, 1], grid_res),
                       np.linspace(bounds[1, 0], bounds[1, 1], grid_res))

    # Construct set of coordinates
    XX = X.flatten()
    YY = Y.flatten()
    poly_arr = np.c_[np.ones(XX.shape[0]), XX, YY, XX * YY, XX * XX, YY * YY]
    Z = np.dot(poly_arr, C)
    para_coords = gl.affine(A, np.c_[XX, YY, Z])

    grid = np.arange(0, XX.shape[0]).reshape((grid_res, grid_res))
    trig_list = np.zeros((6 * (grid_res - 1)**2), dtype=np.int64)

    for j in np.arange(0, grid_res - 1):
        for i in np.arange(0, grid_res - 1):

            # Upper triangle
            trig_list[6 * i + 6 * (grid_res - 1) * j] = grid[(j, i)]
            trig_list[6 * i + 1 + 6 * (grid_res - 1) * j] = grid[(j, i + 1)]
            trig_list[6 * i + 2 + 6 * (grid_res - 1) * j] = grid[(j + 1,
                                                                  i + 1)]

            # Lower Triangle
            trig_list[6 * i + 3 + 6 * (grid_res - 1) * j] = grid[(j + 1, i)]
            trig_list[6 * i + 4 + 6 * (grid_res - 1) * j] = grid[(j, i)]
            trig_list[6 * i + 5 + 6 * (grid_res - 1) * j] = grid[(j + 1,
                                                                  i + 1)]

    return SurfaceMesh(grid.flatten(), para_coords, trig_list.reshape((-1, 3)))


def main():

    p = argparse.ArgumentParser(description="Generate QC pages based "
                                "on a realistic head model .msh "
                                "file and a centroid for "
                                "projection")

    p.add_argument('msh', type=str, help="Path to .msh realistic head model")
    p.add_argument('centroid',
                   type=str,
                   help="Path to .txt centroid coordinates")
    p.add_argument('dscalar',
                   type=str,
                   help="Path to .dscalar file for brain texture")
    p.add_argument('out_html',
                   type=str,
                   help="Path to HTML output file for QC")

    args = p.parse_args()

    msh = args.msh
    centroid = args.centroid
    texture = np.load(args.dscalar)
    out = args.out_html

    # Load in mesh surfaces
    logging.info(f"Loading in brain surfaces from {msh}...")
    brain = decompose_gmsh(msh, GM_ENTITIES)
    domain = QuadraticDomain(HeadModel(msh), np.genfromtxt(centroid))

    logging.info("Discretizing parameteric mesh...")
    C, iR, bounds = domain.C, domain.iR, domain.bounds
    param_mesh = construct_parameteric_mesh(C, iR, bounds)
    gen_mshplot_html(brain, texture, param_mesh, out)


if __name__ == '__main__':
    main()
