#!/usr/bin/env python


'''
Generate a parameteric surface for a given set of points in RAS

Usage:
    parameterize_surface_patch <patch> <norm> <out_prefix>

Arguments:
    <patch>                             Coordinates for a surface patch
    <norm>                              Mean normal vector of surface patch
    <out_prefix>                        Prefix to output files
'''

from docopt import docopt
from fieldopt import geolib as gl
import numpy as np

def main():

    args = docopt(__doc__)

    patch_file      =   args['<patch>']
    norm_file       =   args['<norm>']
    out_prefix      =   args['<out_prefix>']

    #Load input files
    patch = np.fromfile(patch_file).reshape((-1,3))
    vnorm = np.fromfile(norm_file).reshape((-1,3))
    vnorm = vnorm/np.linalg.norm(vnorm)

    #Rotate patch into standard space
    z = np.array([0,0,1], dtype=np.float64)
    R = gl.rotate_vec2vec(vnorm,z)
    r_patch = np.matmul(R,patch.T).T
    
    #Fit quadratic model to Z values
    C = gl.quad_fit(r_patch[:,:2], r_patch[:,2])

    #Make inverse rotation
    inv_R = np.linalg.pinv(R)

    #Minmax bounds
    minarr = np.min(r_patch,axis=0)
    maxarr = np.max(r_patch,axis=0)
    bounds = np.c_[minarr.T,maxarr.T]

    #Save files
    C.tofile(out_prefix + "_C")
    inv_R.tofile(out_prefix + "_R")
    bounds.tofile(out_prefix + "_bounds")

if __name__ == '__main__':
    main()




