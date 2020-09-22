#!/usr/bin/env python

import argparse
import os
import numpy as np
from fieldopt.objective import FieldFunc
from simnibs.msh import mesh_io


def main():

    p = argparse.ArgumentParser(description="Run a SimNIBS simulation "
                                "at the coordinates specified at a given "
                                "orientation")
    p.add_argument("mesh",
                   type=str,
                   help="Realistic head model "
                   ".msh file")
    p.add_argument("orientation",
                   type=str,
                   help="Input parameters to evaluate "
                   "objective function on")
    p.add_argument("centroid",
                   type=str,
                   help="Coordinate .txt file pointing "
                   "to centroid for parameteric mesh seeding")
    p.add_argument("weights",
                   type=str,
                   help="Weight function used to "
                   "evaluate the objective function")
    p.add_argument("coil",
                   type=str,
                   help="Coil to use for running a "
                   "simulation, .ccd physical model or .nii.gz dA/dt file")
    p.add_argument("out_fields",
                   type=str,
                   help="Output gmsh simulation .msh file")
    p.add_argument("out_coil",
                   type=str,
                   help="Output coil position .geo file")
    p.add_argument("out_coords",
                   type=str,
                   help="Output transformed coil coordinate file in RAS")

    args = p.parse_args()

    # Parse arguments
    msh = args.mesh
    orientation = np.genfromtxt(args.orientation)
    centroid = np.genfromtxt(args.centroid)
    wf = np.load(args.weights)
    coil = args.coil
    out_fields = args.out_fields
    out_coil = args.out_coil
    out_coords = args.out_coords

    # Construct the objective function
    fem = FieldFunc(msh,
                    initial_centroid=centroid,
                    tet_weights=wf,
                    coil=coil,
                    field_dir=os.getcwd(),
                    cpus=2)

    # Get position and anterior facing direction of coil
    matsimnibs = fem._transform_input(*orientation)
    a = matsimnibs[:3, 1]

    # Coil adjustment methodology
    # We never want anterior facing coil orientations
    if a[1] < 0:
        orientation[2] = (orientation[2] + 180) % 360
    _, matsimnibs = fem.run_simulation(orientation, out_fields, out_coil)

    # Save matsimnibs matrix
    np.save(out_coords, matsimnibs)

    # Write in weight function
    M = mesh_io.read_msh(out_fields)
    gm = np.where(M.elm.tag1 == 2)
    wf_field = np.zeros_like(M.elmdata[1].value)
    wf_field[gm] = wf
    M.add_element_field(wf_field, 'weightfunction')
    M.write(out_fields)


if __name__ == "__main__":
    main()
