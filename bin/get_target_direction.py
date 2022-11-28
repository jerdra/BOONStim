import argparse
import numpy as np
import nibabel as nib
import fieldopt.geometry.mesh as mesh
from scipy.stats import rankdata

import pyvista as pv
pv.start_xvfb()


def closest_node(coordinate, vertices):
    """
    Find closest node to `coordinate` in the [N,3] vertices array
    """

    diffs = np.linalg.norm(vertices - coordinate, axis=1)
    ind = np.argmin(diffs)
    return diffs[ind], ind


def main():

    p = argparse.ArgumentParser(description="Compute normal direction "
                                "of target coordinate")
    p.add_argument("coordinate", help="Target coordinate", type=str)
    p.add_argument("fs_dir", help="Subject freesurfer directory", type=str)
    p.add_argument("output", help="Output file ending in .npy", type=str)
    p.add_argument("--qc-html", help="Output interactive QC file", type=str)
    p.add_argument("--qc-img", help="Output static QC image", type=str)

    args = p.parse_args()

    l_verts, l_trigs = nib.freesurfer.read_geometry(
        f"{args.fs_dir}/surf/lh.pial")
    l_curv = nib.freesurfer.read_morph_data(f"{args.fs_dir}/surf/lh.curv")

    r_verts, r_trigs = nib.freesurfer.read_geometry(
        f"{args.fs_dir}/surf/rh.pial")
    r_curv = nib.freesurfer.read_morph_data(f"{args.fs_dir}/surf/rh.curv")

    coordinate = np.loadtxt(args.coordinate, delimiter=',')

    l_min, l_ind = closest_node(coordinate, l_verts)
    r_min, r_ind = closest_node(coordinate, r_verts)

    if l_min < r_min:
        ind, verts, trigs, curv = l_ind, l_verts, l_trigs, l_curv
    else:
        ind, verts, trigs, curv = r_ind, r_verts, r_trigs, r_curv

    # Find minimum curvature node to target site
    two_ring = mesh.get_ring(ind, trigs, 2)
    min_curv_ind = np.argmin(np.abs(curv[two_ring]))
    min_curv_vert = two_ring[min_curv_ind]

    # Use minimum curvature node and subset surface for visualization
    proximity_ring = np.sort(mesh.get_ring(min_curv_vert, trigs, 20))
    relevant_triangles = mesh.get_subset_triangles(proximity_ring.astype(int),
                                                   trigs.astype(int))
    proximity_triangles = trigs[np.where(relevant_triangles)]

    # Now crop mesh
    proximity_coords = verts[proximity_ring, :]
    display_inds = np.arange(0, proximity_coords.shape[0])
    display_trigs = rankdata(proximity_triangles, method='dense') - 1
    display_trigs = display_trigs.reshape(proximity_triangles.shape)
    display_curv = curv[proximity_ring]

    # Get two-ring around minimum curvature to compute normal
    _, ind = closest_node(min_curv_vert, proximity_coords)
    normal_ring = mesh.get_ring(ind, display_trigs, 5)
    normal = mesh.get_normals(normal_ring, display_inds, proximity_coords,
                              display_trigs)
    np.save(args.output, -normal)

    # Create display and save
    if args.qc_img or args.qc_html:
        face_arr = np.zeros(
            (display_trigs.shape[0], display_trigs.shape[1] + 1), dtype=int)
        face_arr[:, 0] = 3
        face_arr[:, 1:] = display_trigs
        face_arr = face_arr.flatten()

        surf = pv.PolyData(proximity_coords, face_arr)
        surf.point_data['curvature'] = np.clip(display_curv,
                                               a_min=-0.5,
                                               a_max=0.5)

        opacity_arr = np.ones_like(display_curv) * 0.65
        opacity_arr[normal_ring] = 0
        surf.point_data['opacity'] = opacity_arr

        plotter = pv.Plotter(polygon_smoothing=True, off_screen=True)
        plotter.add_mesh(surf,
                         cmap="RdBu",
                         scalars='curvature',
                         opacity='opacity',
                         use_transparency=True,
                         smooth_shading=True)

        arrow = pv.Arrow(start=coordinate, direction=normal, scale=15)
        sphere = pv.Sphere(radius=1, center=coordinate)
        plotter.add_mesh(arrow, color="black")
        plotter.add_mesh(sphere, color="white")
        plotter.enable_anti_aliasing()
        plotter.window_size = [1000, 1000]

        plotter.camera_position = coordinate + (normal * 20)
        plotter.camera.zoom = 0.50

        if args.qc_img:
            plotter.screenshot(args.qc_img)

        if args.qc_html:
            plotter.export_html(args.qc_html)


if __name__ == '__main__':
    main()
