#!/usr/bin/env python
import argparse
import logging
from collections import namedtuple

import numpy as np
import nibabel as nib
from fieldopt import geolib as gl

from nilearn import image as img
import meshplot as mp

HEAD_ENTITIES = [(2, 5), (2, 1005)]

SurfaceMesh = namedtuple("SurfaceMesh", ["nodes", "coords", "triangles"])
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

    return np.concatenate([surf['left'], surf['right']])


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
        prev_nodes += surf.nodes.max() + 1

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
            return trigs[0].reshape(-1, 3)

    logging.error("Could not properly load Mesh! Check entity tags!")
    raise ValueError


def get_surf_normals(surf, roi):
    '''
    Compute triangle weighted normals for each vertex within an ROI

    Arguments:
        f_cifti:    Path to cifti file
        f_mask:     Path to mask ROI file to constrain normal computation


    Returns:
        roi_inds:   List of ROI vertices for which normals are provided
        norm_array: NX3 array of normals
    '''

    roi_inds = np.where(roi > 0)
    norm_array = np.zeros((len(roi_inds[0]), 3), dtype=np.float)
    all_tags = np.arange(0, surf.coords.shape[0])

    for i, r in enumerate(roi_inds[0]):
        norm_array[i, :] = gl.get_normals([r], all_tags, surf.coords,
                                          surf.triangles)
    return roi_inds[0], norm_array


def get_weighted_centroid(surf, weights):
    '''
    Compute a weighted euclidean centroid given a set of surface
    coordinates and a weighting array

    Arguments:
        surf                SurfaceMesh object
        weights             Array of weights matching the number of vertices
                            in the SurfaceMesh

    Note: Weights will be normalized based on sum to avoid illegitimate
    centroids
    '''

    w = weights / weights.sum()

    centroid = (surf.coords * w).sum(axis=0)
    return centroid


def get_radial_dists(origin_surf, target_surf, roi=None):
    '''
    Given a surface, compute the radial distance projection
    for each of the vertices within the source_roi

    Arguments:
        source:         Surface to generate radial rays from
        target:         Surface to measure distance to
        roi:            Mask surface to constrain regions
                        in which the radial distance will be calculated

    Returns:
        normals:        NX3 array of normals for each of the ROI indices
        ray_dists:      NX1 array of ray distances for each of the ROI indices
    '''

    all_tags = np.arange(0, origin_surf.coords.shape[0])
    if roi is not None:
        roi_inds = np.where(roi > 0)[0]
    else:
        roi_inds = all_tags

    num_coords = len(roi_inds)
    normals = np.zeros_like(origin_surf.coords[roi_inds, :])
    ray_dists = np.zeros((normals.shape[0], 1))

    for i, r in enumerate(roi_inds):
        logging.info(f"Computing radial for vertex {i}/{num_coords}")

        # Generate ray from origin vertex
        n = gl.get_normals(np.array([r]), all_tags, origin_surf.coords,
                           origin_surf.triangles)
        pn = origin_surf.coords[r, :]
        pf = pn + n

        # Compute ray distance
        _, ray_len, _ = gl.ray_interception(pn, pf, target_surf.coords,
                                            target_surf.triangles)

        # Store in output array
        normals[i, :] = n
        ray_dists[i] = ray_len

    return normals, ray_dists


def main():
    parser = argparse.ArgumentParser(description="Given a dscalar image "
                                     "find the head coordinate that "
                                     "approximates a projection of the "
                                     "centroid to the head surface")

    parser.add_argument('mesh', type=str, help="Path to subject .msh model")

    parser.add_argument(
        'dscalar',
        type=str,
        help="Dscalar to find head surface centroid projection of")

    parser.add_argument('left_surf', type=str, help='GIFTI of left hemisphere')

    parser.add_argument('right_surf',
                        type=str,
                        help='GIFTI of right hemisphere')

    parser.add_argument('output_file',
                        type=str,
                        help='Path of text file to output with coordinates')

    args = parser.parse_args()
    f_msh = args.mesh
    f_dscalar = args.dscalar
    f_left = args.left_surf
    f_right = args.right_surf
    output = args.output_file

    logging.info("Loading in surface files...")

    logging.info("gmsh...")
    head = decompose_gmsh(f_msh, HEAD_ENTITIES)

    logging.info("gifti...")
    pial_mesh = compose_surfs(decompose_surf(f_left), decompose_surf(f_right))

    logging.info("dscalar...")
    dscalar = decompose_dscalar(f_dscalar)

    # Compute weighted euclidean centroid
    logging.info("Calculating Euclidean centroid")
    eu_centroid = get_weighted_centroid(pial_mesh, dscalar)

    # Compute normal balanced by strength of weight and ray length
    logging.info("Extracting approximate normal")
    normals, rays = get_radial_dists(pial_mesh, head, dscalar)

    # Punish very far rays significantly
    i_rays = (1/rays) ** 2
    projection_normal = (normals * (i_rays / i_rays.sum())).sum(axis=0)

    # Apply projection normal to eu_centroid
    logging.info("Finding projection of Euclidean centroid to Head")
    pf = eu_centroid + projection_normal
    p_I, _, _ = gl.ray_interception(eu_centroid, pf, head.coords,
                                    head.triangles)

    logging.info(f"Saving results into {output}")
    np.savetxt(output, p_I[:, np.newaxis])


if __name__ == '__main__':
    main()
