#!/usr/bin/env python
import argparse
from fieldopt import geolib as gl
import numpy as np
from nilearn import image as img
import logging

# CONSTS
HEAD_ENTITY = (2, 5)
GM_ENTITY = (2, 2)

logging.basicConfig(format="%(asctime)s [BOONSTIM CORTEX2SCALP]:  %(message)s",
                    datefmt="%Y-%m-%d %I:%M:%S %p",
                    level=logging.INFO)


def decompose_dscalar(dscalar):

    cifti = img.load_img(dscalar)
    data = cifti.get_fdata(dtype=np.float32)
    ax = cifti.header.get_axis(1)
    n_maps = cifti.header.get_axis(0).name.shape[0]

    surf = {}
    for name, indices, brainmodel in ax.iter_structures():

        # Step 1: Get name map
        if not name.startswith("CIFTI_STRUCTURE_CORTEX"):
            continue

        h = name.replace("CIFTI_STRUCTURE_CORTEX_", "").lower()

        # Step 2: Map into structure
        texture = np.zeros((brainmodel.vertex.max() + 1, n_maps),
                           dtype=np.float32)
        texture[brainmodel.vertex, :] = data[:, indices].T
        surf[h] = texture

    return surf


def coil2head(head_coords, f_coil):
    '''
    Project the coil coordinate to the head using minimal euclidean distance
    '''

    c_centre = np.genfromtxt(f_coil)
    head_ind = np.argmin(np.linalg.norm(head_coords - c_centre, axis=1))
    return head_coords[head_ind, :][np.newaxis, :]


def get_min_scalp2cortex(brain_coords, head_coords, f_coil=None):
    '''
    brain_coords: NX3 array containing candidate brain coordinates
    head_coords: MX3 array containing candidate head coordinates
                to check against
    '''

    if f_coil:
        head_coords = coil2head(head_coords, f_coil)

    return get_min_dist(brain_coords, head_coords)


def get_min_dist(x, y, return_coords=False):
    '''
    x: NX3 array
    y: MX3 array

    Find minimial distance between two sets of points

    Returns:
        If return_coords == False:
            minimum distance between specified point sets
        If return_coords == True:
            A tuple of:
                -   minimum distance between specified point sets
                -   coordinate from x that constitutes the minimum
                -   coordinate from y that constitutes the minimum
    '''

    # Calculate set of all pairwise differences
    diffs = x[:, np.newaxis, :] - y
    norms = np.linalg.norm(diffs, axis=2)

    if not return_coords:
        return norms.flatten().min()

    min_ind = np.argmin(norms)

    # Yields row/column combination from x and y
    x_min, y_min = np.unravel_index(min_ind, norms.shape)

    return (norms, x[x_min, :], y[y_min, :])


def apply_mask_to_pial(pial_surf, mask):
    '''
    pial_surf:      Coordinates for pial surface
    mask:           List of mask indices
    '''

    left_mask = np.where(mask['left'] > 0)[0]
    right_mask = np.where(mask['right'] > 0)[0]

    return {
        'left': pial_surf['left'][left_mask, :],
        'right': pial_surf['right'][right_mask, :]
    }


def construct_dist_qc_view(pial_surf, head_surf, coil_centre, mt_roi,
                           geo_file):
    '''
    Function to generate a QC gmsh model

    pial_surf:      Coordinates for pial surface
    head_surf:      Coordinates for head surface
    coil_centre:    Coil centre coordinates
    roi_mask:       Set of coordinates describing MT
    geo_file:       Path to .geo gmsh template file to fill


    NOTE:
        pial_surf, head_surf, coil_centre - must all be in RAS,
        freesurfer T1 native space
    '''

    # Step 1: Project coil to head
    coil_coord = coil2head(head_surf, coil_centre)

    # Step 2: Minimize distance from mt_roi to head
    min_mt, min_mt_head, min_roi = get_min_dist(head_surf,
                                                mt_roi,
                                                return_coords=True)

    # Step 3: Minimize distance from coil to cortex
    min_coil, min_coil_head, min_cortex = get_min_dist(coil_coord,
                                                       pial_surf,
                                                       return_coords=True)

    # configure scalar line
    mt_line = create_scalar_line(min_mt_head, min_roi)
    mt_text = create_text_at_coord((min_mt_head + min_roi) / 2, min_mt)

    coil_line = create_scalar_line(min_coil_head, min_cortex)
    coil_text = create_text_at_coord((min_coil_head + min_cortex) / 2,
                                     min_coil)

    elems = [mt_line, mt_text, coil_line, coil_text]
    view = construct_view("cortex2scalp", elems)

    return view


def merge_mesh(mesh):
    '''
    Given a geo file append mesh merging

    mesh:   Path to gmsh .msh file
    '''

    return f"Merge !{mesh}"


def create_scalar_line(p1, p2):
    '''
    p1: Starting coordinates for scalar line (RAS)
    p2: Ending coordinates for scalar line (RAS)
    '''

    return f'''
    SL({p1[0]},{p1[1]},{p1[2]},{p2[0]},{p2[1]},{p2[2]}){{0,0,0,0,0,0}};
    '''


def create_text_at_coord(p, text):
    '''
    p: Coordinate to insert text
    text: text to insert
    '''

    return f'''
    T3({p[0]},{p[1]},{p[2]},7){{text}};
    '''


def construct_view(name, elems):
    '''
    Given a list of elements to draw, construct a view
    '''

    insert_elems = "\n".join(elems)
    return f'''
    View"{name}"{{
        {insert_elems}
    }};
    '''


def write_geo(cmd, out, opt=None):
    '''
    cmd:    Main set of geo script commands
    out:    Output path
    opt:    Optional arguments for setting up a view
    '''

    write_opt = []
    if opt:
        with open(opt, 'r') as f:
            write_opt = f.readlines()

    with open(out, 'w') as f:
        f.write(cmd)
        f.writelines(write_opt)

    return


def main():

    parser = argparse.ArgumentParser(description="Given a head model and "
                                     "a set of ROI coordinates compute "
                                     "the scalp to cortex distance")

    parser.add_argument('mesh', type=str, help="path to GMSH mesh file")
    parser.add_argument('pial_surface',
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
    parser.add_argument('--qc-geo',
                        type=str,
                        help="Path to geo file to generate a QC model")
    parser.add_argument('--geo-view',
                        type=str,
                        help="Template containing gmsh script instructions "
                        "to use by default"

    args = parser.parse_args()
    f_mesh = args.mesh
    f_coords = args.pial_surface
    f_mask = args.mask
    f_coil = args.coilcentre
    f_out = args.output
    f_qc = args.qc_geo
    f_opt = args.geo_view

    # Step 1: Load in the pial surface
    logging.info("Parsing coordinate CIFTI file...")
    surf_coords = decompose_dscalar(f_coords)

    # Apply mask if available
    if f_mask:
        logging.info("Mask file supplied! Reducing search space...")
        surf_mask = decompose_dscalar(f_mask)
        surf_masked = apply_mask_to_pial(surf_coords, surf_mask)

        if not f_qc:
            coords = np.vstack([surf_masked['left'], surf_masked['right']])
    else:
        coords = np.vstack([surf_coords['left'], surf_coords['right']])

    # Get head model coordinates
    logging.info("Loading gmsh entity for head model")
    _, h_coords, _ = gl.load_gmsh_nodes(f_mesh, HEAD_ENTITY)

    # Compute cortical distance, unless in qc mode
    if f_qc:
        logging.info("Running QC routine")
        masked_coords = np.vstack([surf_masked['left'], surf_masked['right']])
        construct_dist_qc_view(coords, h_coords, f_coil, masked_coords, f_qc)
        geo_out = "\n".join([construct_dist_qc_view, merge_mesh(f_mesh)])

        logging.info(f"Writing geo file to {f_out}")
        write_geo(geo_out, f_out, f_opt)

    else:
        dist = get_min_scalp2cortex(coords, h_coords, f_coil)
        np.save(f_out, dist)


if __name__ == '__main__':
    main()
