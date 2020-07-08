#!/usr/bin/env python
import argparse
from fieldopt import geolib as gl
import numpy as np

HEAD_ENTITY = (2, 5)
GM_ENTITY = (2, 3)


def cortex2scalp(f_mesh, f_roi, ENTITY):
    '''
    Get distance from ROI centroid to closest scalp surface coordinate
    '''

    roi = np.genfromtxt(f_roi)
    h_nodes, h_coords, _ = gl.load_gmsh_nodes(f_mesh, ENTITY)
    return np.linalg.norm(roi - h_coords, axis=1).min()


def scalp2cortex(f_mesh, f_coil, ENTITY):
    '''
    Get distance from Coil centre to closest GM surface coordinate
    '''

    c_centre = np.genfromtxt(f_coil)
    g_nodes, g_coords, _ = gl.load_gmsh_nodes(f_mesh, ENTITY)
    return np.linalg.norm(c_centre - g_coords, axis=1).min()


def main():

    parser = argparse.ArgumentParser(description="Given a head model and "
                                     "a target ROI, compute the distance to "
                                     "the scalp")

    parser.add_argument('mesh', type=str, help="path to GMSH mesh file")

    parser.add_argument('--roi',
                        type=str,
                        help="path to centroid coordinate file")
    parser.add_argument('--coilcentre',
                        type=str,
                        help="path to coil centre coordinate file")
    parser.add_argument('output',
                        type=str,
                        help="output file path containing "
                        "the minimum distance from brain to "
                        "scalp")

    args = parser.parse_args()
    f_mesh = args.mesh
    f_roi = args.roi
    f_out = args.output
    f_coil = args.coilcentre

    if f_roi:
        distance = cortex2scalp(f_mesh, f_roi, HEAD_ENTITY)
    elif f_coil:
        distance = scalp2cortex(f_mesh, f_coil, GM_ENTITY)

    np.save(f_out, distance)


if __name__ == '__main__':
    main()
