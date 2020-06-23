#!/usr/bin/env python
import argparse
from fieldopt import geolib as gl
import numpy as np

HEAD_ENTITY = (2, 5)


def main():

    parser = argparse.ArgumentParser(description="Given a head model and "
                                     "a target ROI, compute the distance to "
                                     "the scalp")

    parser.add_argument('mesh', type=str, help="path to GMSH mesh file")

    parser.add_argument('roi',
                        type=str,
                        help="path to centroid coordinate file")

    parser.add_argument('output',
                        type=str,
                        help="output file path containing "
                        "the minimum distance from brain to "
                        "scalp")

    args = parser.parse_args()
    f_mesh = args.mesh
    f_roi = args.roi
    f_out = args.output

    # Load in ROI coordinate
    roi = np.genfromtxt(f_roi)

    # Load in head
    h_nodes, h_coords, _ = gl.load_gmsh_nodes(f_mesh, HEAD_ENTITY)

    # Find the minimum distance
    scalp2cortex = np.linalg.norm(roi - h_coords, axis=1).min()
    np.save(f_out, scalp2cortex)


if __name__ == '__main__':
    main()
