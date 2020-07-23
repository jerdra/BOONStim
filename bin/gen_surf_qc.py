#!/usr/bin/env python

import argparse
import logging

from nilearn import plotting as plot

import nibabel as nib
import numpy as np

logging.basicConfig(format="%(asctime)s [BOONSTIM GRID]:  %(message)s",
                    datefmt="%Y-%m-%d %I:%M:%S %p",
                    level=logging.INFO)

hemi2vert = {'left': 0, 'right': 1}


def construct_cifti_surf_mesh(l, r):
    cifti_coord = []
    cifti_trigs = []

    prev_max = 0
    for h in [l, r]:
        coords = h.darrays[0].data
        cifti_coord.append(coords)
        cifti_trigs.append(h.darrays[1].data + prev_max)
        prev_max = coords.shape[0]

    return (np.concatenate(cifti_coord,
                           axis=0), np.concatenate(cifti_trigs, axis=0))


def construct_map_from_cifti(cifti):
    vertvals = np.asanyarray(cifti.dataobj)[0]
    brain_models = list(cifti.header.get_index_map(1).brain_models)
    vertices = np.zeros(sum(m.surface_number_of_vertices
                            for m in brain_models),
                        dtype=vertvals.dtype)

    prev_verts = 0
    for b in brain_models:
        inds2fill = np.array(b.vertex_indices) + prev_verts
        vertices[inds2fill] = vertvals[b.index_offset:b.index_offset +
                                       b.index_count]
        prev_verts = b.surface_number_of_vertices

    return vertices


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('dscalar', type=str, help="CIFTI dscalar file to view")
    parser.add_argument("l_surf",
                        type=str,
                        help="Left surface to visualize over")
    parser.add_argument("r_surf",
                        type=str,
                        help="Right surface to visualize over")
    parser.add_argument("output_base",
                        type=str,
                        help="Output basename, will add _view-VIEW to end")
    parser.add_argument("--bg_surf",
                        type=str,
                        help="Optional background surface to modulate "
                        "levels of dscalar map")

    args = parser.parse_args()

    dscalar = args.dscalar
    f_l_surf = args.l_surf
    f_r_surf = args.r_surf
    f_bg = args.bg_surf
    outbase = args.output_base

    # Use background image if available
    if f_bg:
        logging.info(f"Using BG map {f_bg}")
        bg = -construct_map_from_cifti(nib.load(f_bg))
        bg_exists = True
    else:
        logging.info(f"No BG map supplied!")
        bg = None
        bg_exists = False

    l_surf = nib.load(f_l_surf)
    r_surf = nib.load(f_r_surf)

    logging.info("Joining together surface mesh")
    coord, trigs = construct_cifti_surf_mesh(l_surf, r_surf)

    logging.info("Constructing map from dscalar CIFTI")
    vertices = construct_map_from_cifti(nib.load(dscalar))

    # Construct mapping from left/right
    logging.info("Generating anterior view...")
    plot.plot_surf(surf_mesh=[coord, trigs],
                   surf_map=vertices,
                   bg_map=bg,
                   bg_on_data=True if bg_exists else False,
                   darkness=0.8,
                   view='anterior',
                   output_file=f"{outbase}_view-anterior.png")

    logging.info("Generating dorsal view...")
    plot.plot_surf(surf_mesh=[coord, trigs],
                   surf_map=vertices,
                   bg_map=bg,
                   bg_on_data=True if bg_exists else False,
                   darkness=0.8,
                   view='dorsal',
                   output_file=f"{outbase}_view-anterior.png")

    # Hemispheric views
    ind2h = {0: 'left', 1: 'right'}
    prev_vert = 0
    for i, h in enumerate([l_surf, r_surf]):

        # Slice up to length of underlying mesh
        inds = slice(prev_vert, prev_vert + h.darrays[0].data.shape[0])
        h_verts = vertices[inds]
        hemi = ind2h[i]
        h_bg = bg[inds] if bg_exists else None

        logging.info(f"Generating lateral view of {hemi} hemisphere")
        plot.plot_surf(surf_mesh=[h.darrays[0].data, h.darrays[1].data],
                       surf_map=h_verts,
                       bg_map=h_bg,
                       bg_on_data=True if bg_exists else False,
                       darkness=0.8,
                       view='lateral',
                       hemi=hemi,
                       output_file=f"{outbase}_view-lateral_hemi-{hemi}.png")

        logging.info(f"Generating medial view of {hemi} hemisphere")
        plot.plot_surf(surf_mesh=[h.darrays[0].data, h.darrays[1].data],
                       surf_map=h_verts,
                       bg_map=h_bg,
                       bg_on_data=True if bg_exists else False,
                       darkness=0.8,
                       view='medial',
                       hemi=hemi,
                       output_file=f"{outbase}_view-medial_hemi-{hemi}.png")

        prev_vert += h.darrays[0].data.shape[0]


if __name__ == '__main__':
    main()
