import argparse
import json
import numpy as np
from simnibs import opt_struct

from fieldopt.objective import FieldFunc
from fieldopt.geometry.mesh_wrapper import HeadModel


def generate_opt_config(fjson, use_direction):
    """
    Create a SimNIBS TMS optimization config object
    from an input JSON file
    """

    with open(fjson, 'r') as fhandle:
        spec = json.load(fhandle)

    tms_opt = opt_struct.TMSoptimize()

    tms_opt.target = np.array([spec['pos_x'], spec['pos_y'], spec['pos_z']],
                              dtype=float)

    if use_direction:
        tms_opt.target_direction = np.array(
            [spec['dir_x'], spec['dir_y'], spec['dir_z']], dtype=float)

    tms_opt.distance = float(spec['thickness'])

    return tms_opt


def main():
    parser = argparse.ArgumentParser(description="Run ADM optimization")
    parser.add_argument("msh",
                        help="Subject realistic head model .msh file",
                        type=str)
    parser.add_argument("json", help="Optimization settings for msh", type=str)
    parser.add_argument("coil", help=".ccd coil definition file", type=str)
    parser.add_argument("--radius", help="Radius of ROI", type=float)
    parser.add_argument("--sim_result",
                        help="Simulation .msh output",
                        type=str)
    parser.add_argument("--sim_geo", help="Simulation .geo output", type=str)
    parser.add_argument("--ncpus",
                        help="Number of CPUS to use in optimization",
                        type=int,
                        default=4)
    parser.add_argument("orientation",
                        help="Optimal coil orientation matrix output")

    subparsers = parser.add_subparsers(help="Mode of optimization")

    direction_parser = subparsers.add_parser(
        "direction", help="Use E-Field directional optimization")
    direction_parser.set_defaults(direction=True)

    mag_parser = subparsers.add_parser(
        "magnitude", help="Use E-Field magnitude optimization")
    mag_parser.set_defaults(direction=False)

    args = parser.parse_args()
    tms_opt = generate_opt_config(args.json, args.direction)
    tms_opt.fnamehead = args.msh
    tms_opt.fnamecoil = args.coil
    tms_opt.target_size = args.radius
    tms_opt.method = "ADM"
    tms_opt.solver_options = "pardiso"
    tms_opt.angle_resolution = 1
    tms_opt.open_in_gmsh = False

    msn = tms_opt.run(cpus=args.ncpus)
    np.save(args.orientation, msn)

    # Run visualization in FieldOpt, simnibs is difficult..
    model = HeadModel(tms_opt.mesh)
    target_region = tms_opt._get_target_region()

    # Get IDs for GM tetrahedrons
    gm_ids = model.get_tet_ids(2)[0]

    # Shift target regions to be 0 at start of grey matter
    # Subtract by one to account for SimNIBS weirdness
    adjusted_target = target_region - gm_ids.min() - 1

    weightfunc = np.zeros(gm_ids.shape[0])
    weightfunc[adjusted_target] = 1

    func = FieldFunc(model, None, args.coil, tet_weights=weightfunc)
    func.visualize_evaluate(matsimnibs=msn,
                            out_sim=args.sim_result,
                            out_geo=args.sim_geo)


if __name__ == '__main__':
    main()
