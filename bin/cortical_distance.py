#!/usr/bin/env python
import argparse
from fieldopt import geolib as gl
import numpy as np
from nilearn import image as img

# CONSTS
HEAD_ENTITY = (2, 5)
GM_ENTITY = (2, 2)


def decompose_dscalar(dscalar):

    cifti = img.load_img(dscalar)
    data = cifti.get_fdata(np.float32)
    ax = cifti.header.get_axis(1)
    n_maps = cifti.header.get_axis(0).name.shape[0]

    surf = {}
    for name, indices, brainmodel in ax.iter_structures():

        # Step 1: Get name map
        if not name.startswith("CIFTI_STRUCTURE_CORTEX"):
            continue

        h = name.replace("CIFTI_STRUCTURE_CORTEX_", "").lower()

        # Step 2: Map into structure
        texture = np.zeros((n_maps, brainmodel.vertex.max() + 1),
                           dtype=np.float32)
        texture[brainmodel.vertex] = data[:, indices]
        surf[h] = texture

    return surf


def get_min_scalp2cortex(brain_coords, head_coords):
    '''
    brain_coords: NX3 array containing candidate brain coordinates
    head_coords: MX3 array containing candidate head coordinates
                to check against
    '''

    # Get pairwise euclidean distances
    dists = np.linalg.norm(head_coords.reshape[:, np.newaxis, :] -
                           brain_coords,
                           axis=2)

    # Get minimal distance
    return dists.flatten().min()


def main():

    parser = argparse.ArgumentParser(description="Given a head model and "
                                     "a set of ROI coordinates compute "
                                     "the scalp to cortex distance")

    parser.add_argument('mesh', type=str, help="path to GMSH mesh file")
    parser.add_argument('coord_dscalar',
                        type=str,
                        help="path to coordinate dscalar "
                        ",see wb_command --surface-coordinates-to-metric for "
                        "information on how to generate a coordinate dscalar")
    parser.add_argument('--mask',
                        type=str,
                        help="optional mask to speed up computation time")
    parser.add_argument('output',
                        type=str,
                        help="output file path containing "
                        "the minimum distance from brain to "
                        "scalp")
    parser.add_argument('--coilcentre',
                        type=str,
                        help="path to coil centre coordinate file, "
                        "if a coil centre is provided, we first project "
                        "the coordinate to the scalp to calculate effective "
                        "scalp to cortex distance")

    args = parser.parse_args()
    f_mesh = args.mesh
    f_coords = args.coord_dscalar
    f_mask = args.mask
    f_coil = args.coilcentre
    f_out = args.output

    # Step 1: Load in the coordinate dscalar
    surf_coords = decompose_dscalar(f_coords)

    # Apply mask
    if f_mask:
        surf_mask = decompose_dscalar(f_mask)

        # Determine which surfs to maintain
        left_mask = np.where(surf_mask['left'] > 0)
        right_mask = np.where(surf_mask['right'] > 0)

        coords = np.vstack([
            surf_coords['left'][:, left_mask[1]].T,
            surf_coords['right'][:, right_mask[1]].T
        ])
    else:
        coords = np.vstack([surf_coords['left'].T, surf_coords['right'].T])

    # Get head model coordinates
    _, h_coords, _ = gl.load_gmsh_nodes(f_mesh, HEAD_ENTITY)

    # If using coil centre then you must project to the closest head vertex
    if f_coil:
        c_centre = np.genfromtxt(f_coil)
        head_ind = np.argmin(np.linalg.norm(h_coords - c_centre, axis=1))

        # Maintain shape
        h_coords = h_coords[head_ind, :][np.newaxis, :]

    # Now compute the minimal distance from ROI --> mask
    dist = get_min_scalp2cortex(coords, h_coords)

    # Write value into a text file
    np.save(f_out, dist)


if __name__ == '__main__':
    main()
