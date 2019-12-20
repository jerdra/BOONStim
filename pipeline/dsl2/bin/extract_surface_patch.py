#!/usr/bin/env python


'''
Given a .msh file of a tetrahedral head reconstruction and a centroid coordinate in 3D space, generate a patch of surface of scalp that is closest to the centroid coordinate. 

Usage:
    extract_surface_patch.py [options] <msh> <centroid> <out_prefix> 

Arguments:
    <msh>                       Subject .msh file
    <centroid>                  Textfile containing centroid coordinates in RAS space (readable by numpy)
    <out_prefix>                Output binary files prefix containing coordinates for parametric surface
                                

Optional:
    -r, --radius                Radius of expansion around centroid coordinate on scalp [default: 25]
    -d, --distance              Projection distance (for debugging and visualization [default: 1]
'''

import numpy as np
import gmsh
from fieldopt import geolib as gl
from docopt import docopt


# List of possible SURF_HEAD
SURF_HEAD=[(2,1005), (2,5)]
def main():

    arguments = docopt(__doc__)

    msh_file        =   arguments['<msh>']
    centroid_file   =   arguments['<centroid>']
    out             =   arguments['<out_prefix>']
    radius          =   arguments['--radius'] or 25
    distance        =   arguments['--distance'] or 1

    #Try every alternative of SURF_HEAD
    n_tag, n_coord, tri = None, None, None
    for i,s in enumerate(SURF_HEAD):
        try:
            n_tag, n_coord, _ = gl.load_gmsh_nodes(msh_file, s)
            _, _, tri = gl.load_gmsh_elems(msh_file, s)
        except ValueError:
            if i == len(SURF_HEAD):
                raise
            else:
                continue
        else:
            break

    tri = tri[0]
    tri = tri.reshape( (len(tri)//3,-1) )

    #Load in centroid voxel
    centroid = np.genfromtxt(centroid_file)

    #Get minimum euclidean distance
    eudist = np.linalg.norm(n_coord - centroid, axis=1)
    min_ind = np.argmin(eudist)
    
    #Capture nodes within spherical ROI
    head_centroid_to_all = np.linalg.norm(n_coord - n_coord[min_ind], axis=1)
    search_inds = np.where(head_centroid_to_all < radius)

    #Get relevant triangles to vertices
    vert_list = n_tag[search_inds]
    vert_coords = n_coord[search_inds]
    t_arr = gl.get_relevant_triangles(vert_list, tri)
    rel_ind = np.where(t_arr > 0)
    t_rel = tri[rel_ind[0],:]

    #Get triangle to index mapping
    u_val = np.unique(t_rel)
    u_ind = np.arange(0,u_val.shape[0])

    #Map each triangle to its index
    sort_map = {v:i for v,i in zip(u_val,u_ind)}
    map_func = np.vectorize(lambda x: sort_map[x])
    mapped_trigs = map_func(t_rel)
    rel_verts = np.where( np.isin(n_tag,u_val) )
    rel_verts_coords = n_coord[rel_verts,:][0]

    #Compute vertex normals
    norm_arr = gl.get_vert_norms(mapped_trigs, rel_verts_coords)

    #Dilate
    v_norm = np.mean(norm_arr,axis=0)
    dil_coords = vert_coords + distance*v_norm

    #Write dilated vertices and mean normal to file
    np.save(out+"_dilated_coords.npy",dil_coords)
    np.save(out+"_mean_norm.npy",v_norm)

    #Generate param surf
    dil_faces_ind = gl.get_subset_triangles(vert_list, t_rel)
    dil_faces = t_rel[np.where(dil_faces_ind)].flatten(order='C') + vert_list.max()
    dil_faces = list(dil_faces)
    dil_verts = vert_list + vert_list.max()
    dil_coords = dil_coords.flatten()
    gmsh.initialize()
    gmsh.model.add('param_surf')
    tag = gmsh.model.addDiscreteEntity(2, 2001)
    gmsh.model.mesh.setNodes(2, tag, nodeTags=dil_verts, coord=dil_coords)
    gmsh.model.mesh.setElements(2, tag, [2],
                                elementTags=[range(1, len(dil_faces)//3 + 1)],
                                nodeTags=[dil_faces])
    gmsh.write(out+"_param_surf.msh")
    gmsh.finalize()

if __name__ == '__main__':
    main()

