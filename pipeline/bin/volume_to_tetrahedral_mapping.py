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


SURF_HEAD=[(3,1005), (3,5)]
def main():

    args = docopt(__doc__)

    vol_file        =   args['<vol>']
    fem_file        =   args['<FEM>']
    out             =   args['<out>']

    # First load in the image file
    img = nib.load(vol_file)

    #Try every alternative of SURF_HEAD
    for i,s in enumerate(SURF_HEAD):
        try:

        

if __name__ == '__main__':
    main()
