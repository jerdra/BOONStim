#!/usr/bin/env python
import argparse
from fieldopt.geometry import mesh as gl
import numpy as np

HEAD_ENTITY = (2, 5)
GM_ENTITY = (2, 2)


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

    # Load in HEAD model
    _, h_coords, _ = gl.load_gmsh_nodes(f_mesh, HEAD_ENTITY)

    if f_roi:
        roi = np.genfromtxt(f_roi)
        distance = np.linalg.norm(roi - h_coords, axis=1).min()
    elif f_coil:
        # Find closest scalp coordinate from coil
        c_centre = np.genfromtxt(f_coil)
        head_ind = np.argmin(np.linalg.norm(h_coords - c_centre, axis=1))

        # Find scalp to cortex distance here
        _, g_coords, _ = gl.load_gmsh_nodes(f_mesh, GM_ENTITY)
        distance = np.linalg.norm(h_coords[head_ind, :] - g_coords,
                                  axis=1).min()

    np.save(f_out, distance)


if __name__ == '__main__':
    main()
