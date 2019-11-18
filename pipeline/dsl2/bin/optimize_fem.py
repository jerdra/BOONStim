#!/usr/bin/env python

'''
Given a full head model, description of quadratic input domain, and cost weights use Bayesian Optimization from Cornell-MOE to find the optimal input (position and rotation). 

Usage:
    optimize_fem.py [options] <mesh> <weightfile> <quad_const> <bounds> <affine> <coil> <loc_out> <rot_out>

Arguments:
    <mesh>                                  GMSH mesh file (v2 or v4)
    <weightfile>                            1D array consisting of 1 value per node
    <quad_const>                            Quadratic constants for input domain
    <bounds>                                Search bounds (x,y) for quadratic surface
    <affine>                                Affine matrix that takes you from canonical orientation
                                            to surface RAS used in mesh
    <coil>                                  Coil magnetization field (.nii.gz) from SimNIBS
    <loc_out>                               Output position vector in RAS
    <rot_out>                               Output rotation vector to smoothed head surface

Options:
    -c,--cpus N_CPUS                        Number of CPUS to use to calculate objective func
                                            [Default: 8]
    -t,--tmp-dir                            Directory to perform FEM experiments in
                                            [Default: $TMPDIR]
                                            [Fallback: /tmp/]
    -n,--n-iters ITERS                      Maximum number of iterations to perform
                                            [Default: 50]
    -c,--convergence TOL                    Convergence threshold for shifts in input space (decimal number)
                                            [Default: 1e-3]
'''

#Base package loading
import os
import numpy as np
from fieldopt import geolib
from fieldopt.objective import FieldFunc
from docopt import docopt

#Cornell package loading
from moe.optimal_learning.python.cpp_wrappers.domain import TensorProductDomain as cTensorProductDomain
from moe.optimal_learning.python.python_version.domain import TensorProductDomain
from moe.optimal_learning.python.geometry_utils import ClosedInterval
from moe.optimal_learning.python.cpp_wrappers.expected_improvement import ExpectedImprovement
from moe.optimal_learning.python.cpp_wrappers.expected_improvement import multistart_expected_improvement_optimization as meio
from moe.optimal_learning.python.data_containers import HistoricalData, SamplePoint
from moe.optimal_learning.python.cpp_wrappers.log_likelihood_mcmc import GaussianProcessLogLikelihoodMCMC
from moe.optimal_learning.python.default_priors import DefaultPrior
from moe.optimal_learning.python.python_version.optimization import GradientDescentOptimizer, GradientDescentParameters
from moe.optimal_learning.python.cpp_wrappers.optimization import GradientDescentOptimizer as cGDOpt
from moe.optimal_learning.python.cpp_wrappers.optimization import GradientDescentParameters as cGDParams



def gen_sample_from_qei(gp,search_domain,sgd_params,num_samples, num_mc=1e4, lhc_iter=2e4):
        
    qEI = ExpectedImprovement(gaussian_process=gp, num_mc_iterations=int(num_mc))
    optimizer = cGDOpt(search_domain, qEI, sgd_params, int(lhc_iter))
    points_to_sample = meio(optimizer, None, num_samples, use_gpu=False, which_gpu=0,
                            max_num_threads=8)
    qEI.set_current_point(points_to_sample[0])
            
    return points_to_sample, qEI.compute_expected_improvement()

def main():

    args = docopt(__doc__)

    #Parse arguments
    mesh        =   args['<mesh>']
    weights     =   np.load(args['<weightfile>'])
    C           =   np.load(args['<quad_const>'])
    b           =   np.load(args['<bounds>'])
    R           =   np.load(args['<affine>'])
    coil        =   args['<coil>']
    loc_out     =   args['<loc_out>']
    rot_out     =   args['<rot_out>']
    cpus        =   int(args['--cpus']) or 8
    tmpdir      =   args['--tmp-dir'] or os.getenv('TMPDIR') or "/tmp/"
    num_iters   =   args['--n-iters'] or 50
    tol         =   args['--convergence'] or 1e-3


    #Make search domain
    search_domain = TensorProductDomain([
            ClosedInterval(b[0,0],b[0,1]), #X coord on quadratic surface
            ClosedInterval(b[1,0],b[1,1]), #Y coord on quadratic surface
            ClosedInterval(0,180) #Rotational angle
            ])

    c_search_domain = cTensorProductDomain([
            ClosedInterval(b[0,0],b[0,1]), 
            ClosedInterval(b[1,0],b[1,1]),
            ClosedInterval(0,180)
            ])

    #Make objective function
    f = FieldFunc(mesh_file=mesh, quad_surf_consts=C,
                  surf_to_mesh_matrix=R, tet_weights=weights,
                  field_dir=tmpdir, coil=coil, cpus=cpus)


    #Generate historical points
    hist_pts = int(cpus * 1.5)
    init_pts = search_domain.generate_uniform_random_points_in_domain(hist_pts)
    observations = -f.evaluate(init_pts)
    hist_data = HistoricalData(dim = 3, num_derivatives= 0)
    hist_data.append_sample_points([SamplePoint(inp,o,0.0) 
                                    for o,inp in 
                                    zip(observations,init_pts)])

    #Set up model specifications
    prior = DefaultPrior(n_dims = 3 + 2, num_noise=1)
    gp_ll = GaussianProcessLogLikelihoodMCMC(historical_data=hist_data,
                                             derivatives=[], prior=prior,
                                             chain_length=1000, burnin_steps=2000,
                                             n_hypers=2**4, noisy=False)
    gp_ll.train()

    #Initialize grad desc params
    sgd_params = cGDParams(num_multistarts=200, max_num_steps=50,
                           max_num_restarts=2, num_steps_averaged=4,
                           gamma=0.7, pre_mult=1.0, max_relative_change=0.5,
                           tolerance=1.0e-10)

    num_samples = int(cpus*1.3)
    best_point_history = []

    #Fixed number that should never be reached
    prev_min = -9999
    for i in np.arange(0,num_iters):
            
        #Optimize qEI and pick samples
        points_to_sample, ei = gen_sample_from_qei(gp_ll.models[0],
                                                   c_search_domain, sgd_params=sgd_params,
                                                   num_samples=num_samples, num_mc=2**10)

        #Collect observations
        sampled_points = -f.evaluate(points_to_sample)
        evidence = [SamplePoint(c,v,0.0) for c,v in zip(points_to_sample, sampled_points)]

        #Update model
        gp_ll.add_sampled_points(evidence)
        gp_ll.train()

        #Pull model and pull values
        gp = gp_ll.models[0]
        min_point = np.argmin(gp._points_sampled_value)
        min_val = np.min(gp._points_sampled_value)
        best_coord = gp.get_historical_data_copy().points_sampled[min_point]

        print('Recommended Points:')
        print(points_to_sample)
        print('Expected Improvement: {}'.format(ei))
        print('Current Best:')
        print('f(x*)=',min_val)
        print('Coord:', best_coord)

        best_point_history.append(min_val)

        #Check convergence criterion
        if np.abs(prev_min - min_val) < tol:
            print('Convergence reached!')
            print('Previous minimum: {}'.format(prev_min))
            print('Current minimum: {}'.format(min_val))
            print('Tolerance: {}'.format(tol))
            break

    #Once sampling is done take the best point and transform it back into native space
    preaff_loc = geolib.map_param_2_surf(best_coord[0],best_coord[1],C)
    preaff_rot,_ = geolib.map_rot_2_surf(best_coord[0],best_coord[1],best_coord[2],C)
    loc = np.matmul(R,preaff_loc)
    rot = np.matmul(R,preaff_rot)
    np.savetxt(loc_out,loc)
    np.savetxt(rot_out,rot)

if __name__ == '__main__':
    main()
