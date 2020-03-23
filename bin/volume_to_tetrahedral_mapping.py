#!/usr/bin/env python


'''
Given a ribbon constrained volume image with weights per voxel. Project the data into tetrahedral space using the tetrahedral projection algorithm.

Usage:
    volume_to_tetrehedral_mapping <vol> <FEM> <out>

Arguments:
    <vol>                               Ribbon constrained volume file (.nii.gz)
    <FEM>                               GMSH .msh file
    <out>                               Output file (weight per tetrahedron array)

'''

import gmsh
import numpy as np
import nibabel as nib
from docopt import docopt
from fieldopt import geolib as gl
from fieldopt import tetrapro as tp


SURF_HEAD=[(3,1002), (3,2)]


def guess_entity(msh, dim, tag):
    '''
    Use last digit method of figuring out what the entity ID is
    '''

    tag = str(tag)

    gmsh.initialize()
    gmsh.open(msh)
    ent_list = gmsh.model.getEntities()

    subset = [k for k in ent_list if k[0] == dim]
    entity = [k for k in subset if str(k[1])[-1] == tag]

    gmsh.clear()

    return entity[0]

def main():

    args = docopt(__doc__)

    vol_file        =   args['<vol>']
    fem_file        =   args['<FEM>']
    out             =   args['<out>']

    # First load in the image file
    img = nib.load(vol_file)

    tet_entity = guess_entity(fem_file, 3, 2)
    gm_entity = guess_entity(fem_file, 2, 2)
    wm_entity = guess_entity(fem_file, 2, 1)

    tn_tag, tn_coord, _ = gl.load_gmsh_nodes(fem_file,tet_entity)
    te_tag, _, te_param = gl.load_gmsh_elems(fem_file,tet_entity)
    gmn_tag, gmn_coord, _ = gl.load_gmsh_nodes(fem_file,gm_entity)
    wmn_tag, wmn_coord, _ = gl.load_gmsh_nodes(fem_file,wm_entity)

    #Pull ribbon data
    ribbon = img.get_data()
    affine = img.affine
    x,y,z = np.where(ribbon > 0)

    #Formulate lists as arrays
    tn_tag = np.array(tn_tag)
    gmn_tag = np.array(gmn_tag)
    wmn_tag = np.array(wmn_tag)

    #Set up inputs
    tn_list = te_param[0].reshape( (-1,4) )

    #First concatenate nodes
    min_t, max_t, len_t = np.min(tn_tag), np.max(tn_tag), np.size(tn_tag)
    min_g, max_g, len_g = np.min(gmn_tag), np.max(gmn_tag), np.size(gmn_tag)
    min_w, max_w, len_w = np.min(wmn_tag), np.max(wmn_tag), np.size(wmn_tag)

    #Wrap up features
    prop_arr = np.array([
            [min_t, max_t, len_t],
            [min_g, max_g, len_g],
            [min_w, max_w, len_w]
        ])

    #Map nodes to contiguous indexing
    node_list = tp.map_nodes(tn_list, prop_arr)

    #Coordinate array matching indexing
    coord_arr = np.concatenate( (tn_coord, gmn_coord, wmn_coord) )
    coords = coord_arr.reshape((-1,3))

    #Determine whether the file is parcellation-based (integer encoded) or weight-based
    float_encoding = np.any(np.mod(ribbon,1))

    #Project
    if not float_encoding:
        n_out_arr = tp.tetrahedral_parcel_projection(node_list, coord_arr, ribbon, affine)
    else:
        n_out_arr = tp.tetrahedral_weight_projection(node_list, coord_arr, ribbon, affine)
    np.save(out,n_out_arr)


if __name__ == '__main__':
    main()
