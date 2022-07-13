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
from fieldopt.geometry import geometry as gl
import numpy as np

def main():

    args = docopt(__doc__)

    patch_file      =   args['<patch>']
    norm_file       =   args['<norm>']
    out_prefix      =   args['<out_prefix>']

    #Load input files
    patch = np.load(patch_file)
    vnorm = np.load(norm_file)
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
    np.save(out_prefix + "_C",C)
    np.save(out_prefix + "_R",inv_R)
    np.save(out_prefix + "_bounds",bounds)

if __name__ == '__main__':
    main()




