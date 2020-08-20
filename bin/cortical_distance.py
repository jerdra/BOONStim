#!/usr/bin/env python
import argparse
import logging
from collections import namedtuple

import numpy as np
import nibabel as nib
from fieldopt import geolib as gl

from nilearn import image as img
import meshplot as mp

# CONSTS
HEAD_ENTITIES = [(2, 5), (2, 1005)]
GM_ENTITIES = [(2, 2), (2, 1002)]

SurfaceMesh = namedtuple("SurfaceMesh", ["nodes", "coords", "triangles"])
DistResult = namedtuple("DistResult", ["source", "target", "distance"])

logging.basicConfig(format="%(asctime)s [BOONSTIM CORTEX2SCALP]:  %(message)s",
                    datefmt="%Y-%m-%d %I:%M:%S %p",
                    level=logging.INFO)

mp.website()


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
    return get_min_dist(head_coords, brain_coords)


def get_min_dist(x, y):
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
    min_ind = np.argmin(norms)

    # Yields row/column combination from x and y
    x_min, y_min = np.unravel_index(min_ind, norms.shape)

    return DistResult(x[x_min, :], y[y_min, :], norms.flatten().min())


def construct_dist_qcview(s2c, c2s):
    '''
    Function to generate a QC gmsh model

    s2c:            Scalp2Cortex DistResult
    c2s:            Cortex2Scalp DistResult
    geo_file:       Path to .geo gmsh template file to fill


    NOTE:
        pial_surf, head_surf, coil_centre - must all be in RAS,
        freesurfer T1 native space
    '''

    mt_text = create_text_at_coord((c2s.source + c2s.target) / 2, c2s.distance)
    mt_line = create_scalar_line(c2s.source, c2s.target)

    coil_line = create_scalar_line(s2c.source, s2c.target)
    coil_text = create_text_at_coord((s2c.source + s2c.target) / 2,
                                     s2c.distance)

    elems = [mt_line, mt_text, coil_line, coil_text]
    view = construct_view("cortex2scalp", elems)

    return view


def gen_mshplot_html(mesh, out_html, *distresults):
    '''
    Function to generate a QC interactive HTML page rendering

    mesh:           SurfaceMesh object for mesh to render
    *distresults:   DistResult objects to render as lines on top of mesh
    '''

    # Render brain model (will need sulcal information)? think about it
    p = mp.plot(mesh.coords, mesh.triangles)

    source_lines = np.array([d.source for d in distresults])
    target_lines = np.array([d.target for d in distresults])
    p.add_lines(source_lines,
                target_lines,
                shading={
                    "line_width": 10,
                    "line_color": "black"
                })
    p.save(out_html)


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

    to_write = '\n'.join(cmd)
    with open(out, 'w') as f:
        f.write(to_write)
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
            return trigs[0].reshape(-1, 3)

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
    if roi is not None:
        roi_inds = np.where(roi > 0)[0]
    else:
        roi_inds = all_tags

    best_ray = 9999
    best_target = None
    best_source = None
    for i, r in enumerate(roi_inds):

        logging.info(f"Computing radial for vertex {i}/{len(roi_inds)}...")

        # Compute origin surface normal
        n = gl.get_normals(np.array([r]), all_tags, origin_surf.coords,
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

    if best_target is None:
        logging.error("Could not succesfully find cortex to scalp distance!")
        raise ValueError

    return DistResult(best_source, best_target, best_ray)


def write_file(contents, out):
    '''
    Helper func to write contents to out
    '''

    with open(out, 'w') as f:
        f.write(contents)
    return


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
    parser.add_argument('--gmsh-qc',
                        type=str,
                        help="Path to geo file to generate a QC model")
    parser.add_argument('--html-qc',
                        type=str,
                        help="Path to html file to generate an interactive "
                        "qc model")

    args = parser.parse_args()
    f_mesh = args.mesh
    f_lsurf = args.left_surface
    f_rsurf = args.right_surface
    f_roi = args.roi
    f_coil = args.coilcentre
    f_out = args.output
    f_gmsh = args.gmsh_qc
    f_html = args.html_qc

    # TODO: Encode this into argparse
    output_qc = False
    if f_gmsh or f_html:
        output_qc = True
        if not (f_coil and f_roi):
            logging.error(
                "QC mode requires --coilcentre and --roi to be specified!")
            raise ValueError
        scalp2cortex = True
        cortex2scalp = True
    elif f_coil:
        logging.info("Coil centre supplied! "
                     "Calculating scalp under coil to cortex disatance")
        scalp2cortex = True
        cortex2scalp = False
    elif f_roi:
        logging.info("Calculating cortex to scalp distance")
        scalp2cortex = False
        cortex2scalp = True
    else:
        logging.error("Missing critical arguments!")
        logging.error("If cortex2scalp: --roi is needed")
        logging.error("If scalp2cortex: --coilcentre is needed")
        logging.error("If --gmsh-qc: --coilcentre and --roi is needed")
        logging.error("If --html-qc: --coilcentre and --roi is needed")
        raise ValueError

    logging.info("Parsing GIFTI Surface Meshes...")
    pial_mesh = compose_surfs(decompose_surf(f_lsurf), decompose_surf(f_rsurf))

    logging.info("Loading in head model...")
    head = decompose_gmsh(f_mesh, HEAD_ENTITIES)

    if scalp2cortex:
        logging.info("Computing scalp2cortex distance")
        s2c_result = get_min_scalp2cortex(pial_mesh.coords, head.coords,
                                          f_coil)
        if not output_qc:
            write_file(str(s2c_result.distance), f_out)
            return

    if cortex2scalp:
        logging.info("Computing cortex2scalp distance")

        logging.info("Loading ROI file...")
        roi = decompose_dscalar(f_roi)

        logging.info("Computing ray intersection of cortical surface "
                     "to head mesh...")
        c2s_result = get_min_radial_distance(pial_mesh, head, roi=roi)

        if not output_qc:
            write_file(str(c2s_result.distance), f_out)
            return

    if f_gmsh:
        logging.info(f"Writing QC Geo file to {f_out}...")
        geo_contents = construct_dist_qcview(s2c_result, c2s_result)
        write_geo([merge_mesh(f_mesh)] + [geo_contents], f_out, f_gmsh)

    if f_html:
        brain = decompose_gmsh(f_mesh, GM_ENTITIES)
        logging.info(f"Writing QC HTML file to {f_html}...")
        gen_mshplot_html(brain, f_html, s2c_result, c2s_result)


if __name__ == '__main__':
    main()
