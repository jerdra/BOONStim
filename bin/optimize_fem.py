#!/usr/bin/env python
'''
Given a full head model, description of quadratic input domain,
and cost weights use Bayesian Optimization from Cornell-MOE to find the optimal
input (position and rotation).

Usage:
    optimize_fem.py [options] <mesh> <weightfile> <init_centroid> <coil> <output_file>

Arguments:
    <mesh>                                  GMSH mesh file (v2 or v4)
    <weightfile>                            1D array consisting of 1 value per
                                            node
    <init_centroid>                         Initial centroid to start sampling
    <coil>                                  Coil magnetization field (.nii.gz)
                                            from SimNIBS
    <output_file>                           Output basename

Options:
    -c,--cpus N_CPUS                        Number of CPUS to use to calculate
                                            objective func
                                            [Default: 1]
    -t,--tmp-dir TMPDIR                     Directory to perform FEM
                                            experiments in
                                            [Default: $TMPDIR]
                                            [Fallback: /tmp/]
    -n, --n-iters ITERS                     Maximum number of iterations
                                            to perform
                                            [Default: 50]
    -m, --min-var-samps                     Minimum number of samples required
                                            to check convergence
                                            [Default: 10]
    -c,--convergence TOL                    Convergence threshold for shifts
                                            in input space (decimal number)
                                            [Default: 1e-3]
    -h,--history FILE                       Save best point history
                                            [Default: disabled]
    -s,--skip-convergence                   Do not use convergence criterion,
                                            instead use n-iters only
                                            [Default: disabled]
'''

# Base package loading

import os
import logging
from collections import deque
from docopt import docopt

import numpy as np
from fieldopt.objective import FieldFunc

# Cornell package loading

from moe.optimal_learning.python.cpp_wrappers.domain import TensorProductDomain as cTensorProductDomain
from moe.optimal_learning.python.python_version.domain import TensorProductDomain
from moe.optimal_learning.python.geometry_utils import ClosedInterval
from moe.optimal_learning.python.cpp_wrappers.expected_improvement import ExpectedImprovement
from moe.optimal_learning.python.cpp_wrappers.expected_improvement import multistart_expected_improvement_optimization as meio
from moe.optimal_learning.python.data_containers import HistoricalData, SamplePoint
from moe.optimal_learning.python.cpp_wrappers.log_likelihood_mcmc import GaussianProcessLogLikelihoodMCMC
from moe.optimal_learning.python.default_priors import DefaultPrior
from moe.optimal_learning.python.python_version.optimization import GradientDescentOptimizer
from moe.optimal_learning.python.cpp_wrappers.optimization import GradientDescentOptimizer as cGDOpt
from moe.optimal_learning.python.cpp_wrappers.optimization import GradientDescentParameters as cGDParams
from moe.optimal_learning.python.base_prior import TophatPrior, NormalPrior

logging.basicConfig(format="%(asctime)s [BOONSTIM GRID]:  %(message)s",
                    datefmt="%Y-%m-%d %I:%M:%S %p")


def gen_sample_from_qei(gp,
                        search_domain,
                        sgd_params,
                        num_samples,
                        num_mc=1e4,
                        lhc_iter=2e4):

    qEI = ExpectedImprovement(gaussian_process=gp,
                              num_mc_iterations=int(num_mc))
    optimizer = cGDOpt(search_domain, qEI, sgd_params, int(lhc_iter))
    points_to_sample = meio(optimizer,
                            None,
                            num_samples,
                            use_gpu=False,
                            which_gpu=0,
                            max_num_threads=8)
    qEI.set_current_point(points_to_sample[0])

    return points_to_sample, qEI.compute_expected_improvement()


def main():

    args = docopt(__doc__)

    # Parse arguments
    mesh = args['<mesh>']
    weights = np.load(args['<weightfile>'])
    init_centroid = np.genfromtxt(args['<init_centroid>'])
    coil = args['<coil>']
    output_file = args['<output_file>']
    cpus = int(args['--cpus']) or 8
    tmpdir = args['--tmp-dir'] or os.getenv('TMPDIR') or "/tmp/"
    num_iters = int(args['--n-iters']) or 50
    min_samps = int(args['--min-var-samps']) or 10
    tol = float(args['--convergence']) or 0.001
    history = args['--history']
    skip_convergence = args['--skip-convergence']

    logging.info('Using {} cpus'.format(cpus))

    f = FieldFunc(mesh_file=mesh,
                  initial_centroid=init_centroid,
                  tet_weights=weights,
                  coil=coil,
                  field_dir=tmpdir,
                  cpus=cpus)

    # Make search domain
    search_domain = TensorProductDomain([
        ClosedInterval(f.bounds[0, 0], f.bounds[0, 1]),
        ClosedInterval(f.bounds[1, 0], f.bounds[1, 1]),
        ClosedInterval(0, 180)
    ])

    c_search_domain = cTensorProductDomain([
        ClosedInterval(f.bounds[0, 0], f.bounds[0, 1]),
        ClosedInterval(f.bounds[1, 0], f.bounds[1, 1]),
        ClosedInterval(0, 180)
    ])

    # Generate historical points
    prior = DefaultPrior(n_dims=3 + 2, num_noise=1)
    prior.tophat = TophatPrior(-2, 5)
    prior.ln_prior = NormalPrior(12.5, 1.6)
    hist_pts = cpus
    i = 0
    init_pts = search_domain.generate_uniform_random_points_in_domain(hist_pts)
    observations = -f.evaluate(init_pts)
    hist_data = HistoricalData(dim=3, num_derivatives=0)
    hist_data.append_sample_points(
        [SamplePoint(inp, o, 0.0) for o, inp in zip(observations, init_pts)])

    # Train GP model
    gp_ll = GaussianProcessLogLikelihoodMCMC(historical_data=hist_data,
                                             derivatives=[],
                                             prior=prior,
                                             chain_length=1000,
                                             burnin_steps=2000,
                                             n_hypers=2**4,
                                             noisy=False)
    gp_ll.train()

    # Initialize grad desc params
    sgd_params = cGDParams(num_multistarts=200,
                           max_num_steps=50,
                           max_num_restarts=5,
                           num_steps_averaged=4,
                           gamma=0.7,
                           pre_mult=1.0,
                           max_relative_change=0.5,
                           tolerance=1.0e-10)

    num_samples = int(cpus * 1.3)
    best_point_history = []

    # Sum of errors buffer
    var_buffer = deque(maxlen=min_samps)
    for i in np.arange(0, num_iters):

        # Optimize qEI and pick samples
        points_to_sample, ei = gen_sample_from_qei(gp_ll.models[0],
                                                   c_search_domain,
                                                   sgd_params=sgd_params,
                                                   num_samples=num_samples,
                                                   num_mc=2**10)

        # Collect observations
        sampled_points = -f.evaluate(points_to_sample)
        evidence = [
            SamplePoint(c, v, 0.0)
            for c, v in zip(points_to_sample, sampled_points)
        ]

        # Update model
        gp_ll.add_sampled_points(evidence)
        gp_ll.train()

        # Pull model and pull values
        gp = gp_ll.models[0]
        min_point = np.argmin(gp._points_sampled_value)
        min_val = np.min(gp._points_sampled_value)
        best_coord = gp.get_historical_data_copy().points_sampled[min_point]

        logging.info('Iteration {} of {}'.format(i, num_iters))
        logging.info('Recommended Points:')
        logging.info(points_to_sample)
        logging.info('Expected Improvement: {}'.format(ei))
        logging.info('Current Best:')
        logging.info(f'f(x*)= {min_val}')
        logging.info(f'Coord: {best_coord}')
        best_point_history.append(str(min_val))

        if history:
            with open(history, 'w') as buf:
                buf.write('\n'.join(best_point_history))

        # Convergence check
        if (len(var_buffer) == var_buffer.maxlen) and not skip_convergence:
            deviation = sum([abs(x - min_val) for x in var_buffer])
            if deviation < tol:
                logging.info('Convergence reached!')
                logging.info('Deviation: {}'.format(deviation))
                logging.info('History length: {}'.format(var_buffer.maxlen))
                logging.info('Tolerance: {}'.format(tol))
                break

        var_buffer.append(min_val)

    # Save position and orientation matrix
    np.savetxt(output_file, best_coord)


if __name__ == '__main__':
    main()
