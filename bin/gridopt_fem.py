#!/usr/bin/env python

import os
import argparse
import json
import logging

import numpy as np
from sklearn.utils.extmath import cartesian

from fieldopt.objective import FieldFunc

logging.basicConfig(format="%(asctime)s [BOONSTIM GRID]:  %(message)s",
                    datefmt="%Y-%m-%d %I:%M:%S %p")


def main():

    parser = argparse.ArgumentParser(
        description="Run grid optimization on a single subject")
    parser.add_argument('msh',
                        type=str,
                        help="Subject Gmsh .msh realistic head model")
    parser.add_argument('weights',
                        type=str,
                        help=".npy binary containing a weight for each "
                        "tetrahedron")
    parser.add_argument('centroid',
                        type=str,
                        help="Coordinates in T1w space for a centroid "
                        "to the weight function to optimize over")
    parser.add_argument('coil', type=str, help="Path to SimNIBS coil file")
    parser.add_argument('output_file',
                        type=str,
                        help="Output file storing optimal coordinates")
    parser.add_argument('locdim',
                        type=int,
                        help="Number of points to evaluate along each "
                        "spatial dimension")
    parser.add_argument('rotdim',
                        type=int,
                        help="Number of points to evaluate along each "
                        "rotational dimension")
    parser.add_argument('--history',
                        type=str,
                        help="Output file to store history of scores"
                        " into for convergence/visualization")
    parser.add_argument('--workdir',
                        type=str,
                        help="Working directory to run simulations in")
    parser.add_argument('--ncpus',
                        type=int,
                        help="Number of threads to use for each batch "
                        "of simulations. Default = 8")
    parser.add_argument('--batchsize',
                        type=int,
                        help="Number of simulations to run simultaneously, "
                        "will default to half the number of cpus if not "
                        "specified.")
    parser.add_argument('--options',
                        type=str,
                        help="ADVANCED: Modify defaults for FEM evaluation "
                        "function.")

    args = parser.parse_args()
    msh = args.msh
    wf = np.load(args.weights)
    centroid = np.genfromtxt(args.centroid)
    coil = args.coil
    ncpus = args.ncpus or 8
    batch_size = args.batchsize or (ncpus // 2 - 1)
    history = args.history
    workdir = args.workdir or "/tmp/"
    loc_dim = args.locdim
    rot_dim = args.rotdim
    output_file = args.output_file
    options = args.options

    if options:
        with open(options, 'r') as f:
            opts = json.load(f)

    # Construct objective function object
    logging.info(f"Using {ncpus} cpus")
    femfunc = FieldFunc(mesh_file=msh,
                        initial_centroid=centroid,
                        tet_weights=wf,
                        coil=coil,
                        field_dir=workdir,
                        cpus=ncpus,
                        **opts)

    # Set up grid for evaluation
    x_in = np.linspace(femfunc.bounds[0, 0], femfunc.bounds[0, 1], loc_dim)
    y_in = np.linspace(femfunc.bounds[1, 0], femfunc.bounds[1, 1], loc_dim)
    rot_in = np.linspace(0, 180, rot_dim)
    input_array = cartesian([x_in, y_in, rot_in])

    score_list = []
    if history and os.path.exists(history):

        # Clip evaluation array based on number previously evaluated
        logging.info(f"{history} already exists! Skipping previous runs")
        with open(history, 'r') as f:
            num_previous = sum(1 for _ in f) - 1
        input_array = input_array[num_previous:, :]
        score_list.append(np.genfromtxt(history, skip_header=1, delimiter=','))

    elif not os.path.exists(history):
        header = "x,y,r,score\n"
        with open(history, 'w') as f:
            f.write(header)

    # Split input array into chunks
    divisions = np.arange(batch_size, input_array.shape[0], batch_size)
    input_arrays = np.split(input_array, divisions)

    logging.info(f"Running {len(input_arrays)} iterations...")
    logging.info(f"Evaluating {batch_size} simulations per iteration...")

    for i, a in enumerate(input_arrays):

        logging.info(f"Iteration number {i} of {len(input_arrays)}")
        scores = femfunc.evaluate(a)

        score_arr = np.c_[a, scores]
        score_list.append(score_arr)
        if history:
            with open(history, "a") as h_file:
                logging.debug(f"Writing to {history}...")
                np.savetxt(h_file, np.c_[a, scores], delimiter=',')

        logging.info(f"Completed iteration {i}")

    logging.info("Finished computing evaluation grid!")

    # Find the best value
    all_scores = np.vstack(score_list)
    best_row = np.argmax(all_scores[:, -1])
    best_input = all_scores[best_row, :-1]
    np.savetxt(output_file, best_input)
    logging.info(f"Saved optimal coordinates to {output_file}")


if __name__ == '__main__':
    main()
