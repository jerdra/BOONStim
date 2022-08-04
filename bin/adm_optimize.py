import argparse
import numpy as np
from simnibs import opt_struct

from fieldopt.objective import FieldFunc
from fieldopt.geometry.mesh_wrapper import HeadModel


def main():
    p = argparse.ArgumentParser(description="Run ADM optimization")
    p.add_argument("msh",
                   help="Subject realistic head model .msh file",
                   type=str)
    p.add_argument("coordinates",
                   help="Text file of RAS coordinates for target",
                   type=str)
    p.add_argument("distance",
                   help="Distance from scalp to use in optimization",
                   type=float)
    p.add_argument("coil", help=".ccd coil definition file", type=str)
    p.add_argument("--radius", help="Radius of ROI", type=float)
    p.add_argument("--sim_result", help="Simulation .msh output", type=str)
    p.add_argument("--sim_geo", help="Simulation .geo output", type=str)
    p.add_argument("--ncpus",
                   help="Number of CPUS to use in optimization",
                   type=int,
                   default=4)
    p.add_argument("orientation",
                   help="Optimal coil orientation matrix output")

    subparsers = p.add_subparsers(help="Mode of optimization")

    radial_parser = subparsers.add_parser(
        "radial", help="Use E-Field directional optimization")
    radial_parser.add_argument(
        "direction",
        help="Vector text file containing E-Field direction",
        type=str)

    mag_parser = subparsers.add_parser(
        "magnitude", help="Use E-Field magnitude optimization")
    mag_parser.set_defaults(direction=None)

    args = p.parse_args()

    tms_opt = opt_struct.TMSoptimize()
    tms_opt.fnamehead = args.msh
    tms_opt.fnamecoil = args.coil
    tms_opt.target = np.load(args.coordinates)
    if args.direction is not None:
        tms_opt.target_direction = np.load(args.direction)
    tms_opt.distance = args.distance
    tms_opt.target_size = args.radius
    tms_opt.method = "ADM"
    tms_opt.solver_options = "pardiso"
    tms_opt.angle_resolution = 1
    tms_opt.open_in_gmsh = False

    msn = tms_opt.run(cpus=args.ncpus)
    msn.save(args.orientation)

    # Run visualization in FieldOpt, simnibs is difficult..
    model = HeadModel(tms_opt.mesh)
    target_region = tms_opt._get_target_region()

    wf = np.zeros(tms_opt.mesh.elm.nr)
    wf[target_region - 1] = 1

    f = FieldFunc(model, None, args.coil, tet_weights=wf)
    f.visualize_evaluate(matsimnibs=msn,
                         out_sim=args.sim_result,
                         out_geo=args.sim_geo)
