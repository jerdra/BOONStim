import argparse
import numpy as np
import logging
from fieldopt.geometry.domains import QuadraticDomain
from fieldopt.geometry.mesh_wrapper import HeadModel
from fieldopt.objective import FieldFunc

logging.basicConfig(level=logging.INFO)


def get_grid_optimizer(f, args):
    '''
    Retrieve Grid optimizer
    '''
    from fieldopt.optimization import grid

    locations = args.n_locations
    rotations = args.n_rotations
    return grid.get_default_tms_optimizer(f, locations, rotations)


def get_bayes_optimizer(f, args):
    '''
    Retrieve Bayes Optimizer
    '''
    from fieldopt.optimization import bayes_moe

    args = vars(args)
    kwargs = {
        k: args[k]
        for k in ["minimum_samples", "max_iterations"] if args[k] is not None
    }

    return bayes_moe.get_default_tms_optimizer(f, args["nworkers"], **kwargs)


def main():

    parser = argparse.ArgumentParser(description="Run Optimization"
                                     " on a given FEM problem")
    parser.add_argument("mesh", type=str, help="SimNIBS gmsh .msh file")
    parser.add_argument("weightfunction",
                        type=str,
                        help="Weight function 1-D .npy file with 1 node per"
                        " grey matter tetrahedron")
    parser.add_argument(
        "initial_point",
        type=str,
        help="Initial coordinates to construct sampling surface "
        "from")
    parser.add_argument(
        "coil",
        type=str,
        help="Coil to use to run simulations, may be .ccd or .nii.gz")
    parser.add_argument("out_coords",
                        type=str,
                        help="Output path to .txt optimal coordinate file")
    parser.add_argument("--out_msh",
                        type=str,
                        help="Output path of .msh result file")
    parser.add_argument("--out_geo",
                        type=str,
                        help="Output path of .geo coil position file")
    parser.add_argument("--nworkers",
                        type=int,
                        help="Number of workers to spawn per run",
                        default=4)
    parser.add_argument("--history",
                        type=str,
                        help="Path to write history file into")
    parser.add_argument("--ncores",
                        type=int,
                        help="Maximum number of physical cores to use "
                        "when setting up and solving simulations")

    sub_parsers = parser.add_subparsers()
    parser_bayes = sub_parsers.add_parser("bayesian",
                                          help="Use Bayesian Optimization")
    parser_bayes.add_argument(
        "--minimum_samples",
        type=int,
        help="Minimum number of samples to collect before "
        "assessing convergence",
        default=10)
    parser_bayes.add_argument("--max_iterations",
                              type=int,
                              help="Maximum number of iterations to perform")
    parser_bayes.set_defaults(func=get_bayes_optimizer)

    parser_grid = sub_parsers.add_parser("grid", help="Use Grid Optimization")
    parser_grid.add_argument("n_locations",
                             type=int,
                             help="N locations to sample in N x N grid")
    parser_grid.add_argument("n_rotations",
                             type=int,
                             help="N rotations to sample between [0, 180]")
    parser_grid.set_defaults(func=get_grid_optimizer)

    p = parser.parse_args()
    print(p)
    model = HeadModel(p.mesh)
    domain = QuadraticDomain(model, np.genfromtxt(p.initial_point))

    # Set up FEM objective function
    func_kwargs = {
        "head_model": model,
        "sampling_domain": domain,
        "tet_weights": np.load(p.weightfunction),
        "coil": p.coil,
        "nworkers": p.nworkers,
        "nthreads": p.ncores
    }
    f = FieldFunc(**func_kwargs)

    # Optimize and write history
    optimizer = p.func(f, p)
    optimizer.optimize(print_status=True)

    if p.history:
        history = optimizer.get_history()
        np.save(p.history, history)

    # Save optimal coordinates
    i, score = optimizer.current_best
    np.savetxt(p.out_coords, i)

    # Generate best simulation msh and geo
    f.visualize_evaluate(*i, out_sim=p.out_msh, out_geo=p.out_geo)


if __name__ == '__main__':
    main()
