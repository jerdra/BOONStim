#!/usr/bin/env python
import argparse
from fieldopt import geolib as gl
import numpy as np
import nibabel as nib
from nilearn import image as img
import logging
from collections import namedtuple

# CONSTS
HEAD_ENTITIES = [(2, 5), (2, 1005)]
SurfaceMesh = namedtuple("SurfaceMesh", ["ids", "coords", "triangles"])
DistResult = namedtuple("DistResult", ["source", "target", "distance"])

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

        hemi = name.replace("CIFTI_STRUCTURE_CORTEX_", "").lower()

        # Step 2: Map into structure
        texture = np.zeros((brainmodel.vertex.max() + 1, n_maps),
                           dtype=np.float32)
        texture[brainmodel.vertex, :] = data[:, indices].T
        surf[hemi] = texture

    return surf


def decompose_surf(surf):
    '''
    Given a GIFTI surface file return:
        - the index array
        - the coordinate array
        - the triangle array
    '''

    gfti = nib.load(surf)

    coords = gfti.darrays[0].data
    trigs = gfti.darrays[1].data
    node_ids = np.arange(0, coords.shape[0])

    return SurfaceMesh(node_ids, coords, trigs)


def compose_surfs(*surfs):
    '''
    Given a set of SurfaceMesh generate a new
    SurfaceMesh

    Arguments:
        ...SurfaceMesh      Set of mesh to combine

    Output:
        SurfaceMesh
    '''

    prev_nodes = 0
    new_nodes = []
    new_coords = []
    new_trigs = []
    for surf in surfs:
        new_nodes.append(surf.nodes + prev_nodes)
        new_coords.append(surf.coords)
        new_trigs.append(surf.triangles + prev_nodes)
        prev_nodes += surf.nodes.max()

    return SurfaceMesh(np.concatenate(new_nodes), np.vstack(new_coords),
                       np.vstack(new_trigs))


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
    nodes, coords, _ = load_surf_verts(f_mesh, HEAD_ENTITIES)
    trigs = load_surf_trigs(f_mesh, HEAD_ENTITIES)

    return SurfaceMesh(nodes, coords, trigs)


def coil2head(head_coords, f_coil):
    '''
    Project the coil coordinate to the head using minimal euclidean distance
    '''

    c_centre = np.genfromtxt(f_coil)
    head_ind = np.argmin(np.linalg.norm(head_coords - c_centre, axis=1))
    return head_coords[head_ind, :][np.newaxis, :]


def get_min_scalp2cortex(brain_coords, head_coords, f_coil):
    '''
    brain_coords: NX3 array containing candidate brain coordinates
    head_coords: MX3 array containing candidate head coordinates
                to check against
    '''

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

    return DistResult(norms.flatten().min(), x[x_min, :], y[y_min, :])


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

    mt_text = create_text_at_coord((min_mt_head + min_roi) / 2, min_mt)
    mt_line = create_scalar_line(min_mt_head, min_roi)

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

    return f'Merge "{mesh}";\n'


def create_scalar_line(p1, p2):
    '''
    p1: Starting coordinates for scalar line (RAS)
    p2: Ending coordinates for scalar line (RAS)
    '''

    return f'''
    SL({p1[0]},{p1[1]},{p1[2]},{p2[0]},{p2[1]},{p2[2]}){{0,0,0,0,0,0}};
    '''


def create_scalar_point(p):
    '''
    Create sphere at specified coordinate

    p:  Coordinates at which to create sphere
    r:  Radius of sphere to be created
    '''

    return f'''
    SP({p[0]}, {p[1]}, {p[2]}){{0.0}};
    '''


def create_text_at_coord(p, text):
    '''
    p: Coordinate to insert text
    text: text to insert
    '''

    return f'''
    T3({p[0]},{p[1]},{p[2]},7){{"{text:.2f}"}};
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


def load_surf_verts(f_msh, entities):
    '''
    surf:       Path to gmsh MSH file
    entities:   List of entities to attempt of form (dim, tag)
    '''

    for dim, tag in entities:

        try:
            nodes, coords, params = gl.load_gmsh_nodes(f_msh, (dim, tag))
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
            _, _, trigs = gl.load_gmsh_elems(f_msh, (dim, tag))
        except ValueError:
            continue
        else:
            return trigs

    logging.error("Could not properly load Mesh! Check entity tags!")
    raise ValueError


def get_surf_normals(surf, roi):
    '''
    f_cifti:    Path to cifti file
    f_mask:     Path to mask ROI file to constrain normal computation

    Compute triangle weighted normals for each vertex
    '''

    roi_inds = np.where(roi > 0)
    norm_array = np.zeros((len(roi_inds[0]), 3), dtype=np.float)
    all_tags = np.arange(0, surf.coords.shape[0])

    for i, r in enumerate(roi_inds[0]):
        norm_array[i, :] = gl.get_normals([r], all_tags, surf.coords,
                                          surf.triangles)
    return roi_inds[0], norm_array


def get_min_radial_distance(origin_surf, target_surf, roi=None):
    '''
    Given two surfaces, compute the minimum radial distances across
    all vertices between two surface mesh.

    Computation of distances can be restricted using an ROI mask

    Arguments:
        origin_surf             Surface to generate radial rays from
        target_surf             Surface to measure distance to
        roi                     Mask surface to constrain regions
                                in which the radial distance will
                                be calculated
    '''

    all_tags = np.arange(0, origin_surf.coords.shape[0])
    if roi:
        roi_inds = np.where(roi > 0)[0]
    else:
        roi_inds = all_tags

    best_ray = 9999
    best_target = None
    best_source = None
    for r in roi_inds:

        # Compute origin surface normal
        n = gl.get_normals([r], all_tags, origin_surf.coords,
                           origin_surf.triangles)

        # Generate a ray
        pn = origin_surf.coords[r, :]
        pf = pn + n

        # Compute ray distance
        p_I, ray_len, _ = gl.ray_interception(pn, pf, target_surf.coords,
                                              target_surf.triangles)
        if ray_len is None:
            continue
        if ray_len < best_ray:
            best_ray = ray_len
            best_source = pn
            best_target = p_I

    return DistResult(best_source, best_target, best_ray)


def main():

    parser = argparse.ArgumentParser(description="Given a head model and "
                                     "a set of ROI coordinates compute "
                                     "the scalp to cortex distance")

    parser.add_argument('mesh', type=str, help="path to GMSH mesh file")
    parser.add_argument('left_surface',
                        type=str,
                        help="path left surface GIFTI ")
    parser.add_argument('right_surface',
                        type=str,
                        help="path right surface GIFTI ")
    parser.add_argument('--roi',
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

    args = parser.parse_args()
    f_mesh = args.mesh
    f_lsurf = args.left_surface
    f_rsurf = args.right_surface
    f_roi = args.roi
    f_coil = args.coilcentre
    f_out = args.output
    f_qc = args.qc_geo

    # Make sure input args make sense
    if f_qc and f_coil:
        logging.info("Running QC mode!")
        output_qc = True
        scalp2cortex = True
        cortex2scalp = True
    elif f_coil and not f_qc:
        logging.info("Coil centre supplied! "
                     "Calculating scalp under coil to cortex distance")
        output_qc = False
        scalp2cortex = True
        cortex2scalp = False
    elif f_qc and not f_coil:
        logging.error("QC mode requires --coilcentre to be specified!")
        raise ValueError
    else:
        logging.info("Calculating cortex to scalp distance")
        scalp2cortex = False
        cortex2scalp = True

    # Step 1: Load in the pial surfaces
    logging.info("Parsing GIFTI Surface Meshes...")
    pial_mesh = compose_surfs(decompose_surf(f_lsurf), decompose_surf(f_rsurf))

    # Get head model coordinates
    logging.info("Loading in head model...")
    head = decompose_gmsh(f_mesh, HEAD_ENTITIES)

    # Calculate scalp to cortex distance
    if scalp2cortex:
        logger.info("Computing scalp2cortex distance")
        s2c_result = get_min_scalp2cortex(pial_mesh.coords, head.coords,
                                          f_coil)
        print(s2c_result)

    if cortex2scalp:
        logger.info("Computing cortex2scalp distance")
        logger.info("Obtaining cortical surface normals...")
        roi = decompose_dscalar(f_roi)

        logger.info("Computing ray intersection of cortical surface "
                    "to head mesh")

        # Insert routine to compute the ray interesection (geolib func)

    if output_qc:
        pass


if __name__ == '__main__':
    main()
