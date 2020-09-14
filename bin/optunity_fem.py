#!/usr/bin/env python

import os
import argparse
import logging

import numpy as np
from fieldopt.objective import FieldFunc

import optunity

logging.basicConfig(format="%(asctime)s [BOONSTIM GRID]:  %(message)s",
                    datefmt="%Y-%m-%d %I:%M:%S %p")

# How do we evaluate multiple outcomes at once?
def evaluate_objective(x,y,r,f):
    pass

def main():

    parser = argparse.ArgumentParsers(
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
    parser.add_argument('output_file',
                        type=str,
                        help="Output file storing optimal coordinates")
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
    parser.add_argument('--solver',
                        type=int,
                        help="Optunity solver to use, "
                        "defaults to particle swarm",
                        choices=optunity.available_solvers())

    args = parser.parse()
    msh = args.msh
    wf = np.load(args.weights)
    centroid = np.genfromtxt(args.centroid)
    coil = args.coil
    ncpus = args.ncpus or 8
    batch_size = args.batchsize or (ncpus // 2 - 1)
    history = args.history
    workdir = args.workdir or "/tmp/"
    output_file = args.output_file
    solver = args.solver or "particle swarm"

    # Construct objective function object
    f = FieldFunc(mesh_file=msh,
                  initial_centroid=centroid,
                  tet_weights=wf,
                  coil=coil,
                  field_dir=workdir,
                  cpus=ncpus)

    # Set up optunity optimization
    # Can we feed a list of inputs here?
    pars, details, _ = optunity.minimize(f.evaluate,
                                         num_evals=100,
                                         x=[f.bounds[0, 0], f.bounds[0, 1]],
                                         y=[f.bounds[1, 0], f.bounds[1, 1]],
                                         theta=[0, 180],
                                         solver_name=solver)
