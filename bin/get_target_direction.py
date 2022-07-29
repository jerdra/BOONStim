import argparse
import numpy as np
import nibabel as nib
from fieldopt.geometry.mesh import get_ring, get_normals


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

    two_ring = get_ring(ind, trigs, 2)
    min_curv_ind = np.argmin(np.abs(curv[two_ring]))
    min_curv_vert = two_ring[min_curv_ind]

    normal = get_normals(min_curv_vert, range(trigs.shape[0]), verts, trigs)
    normal.save(args.coordinate)


if __name__ == '__main__':
    main()
