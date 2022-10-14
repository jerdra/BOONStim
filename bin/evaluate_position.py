import argparse
import json
import numpy as np
from simnibs import opt_struct

from fieldopt.objective import FieldFunc
from fieldopt.geometry.mesh_wrapper import HeadModel


def get_target_direction(fjson):
    """
    Pull target direction vector from JSON file
    """

    with open(fjson, 'r') as fhandle:
        spec = json.load(fhandle)
    return np.array([spec['dir_x'], spec['dir_y'], spec['dir_z']], dtype=float)


def main():
    parser = argparse.ArgumentParser(description="Run ADM optimization")
    parser.add_argument("msh",
                        help="Subject realistic head model .msh file",
                        type=str)
    parser.add_argument("msn", help="Coil orientation matrix", type=str)
    parser.add_argument("json", help="Optimization settings for msh", type=str)
    parser.add_argument("coil", help=".ccd coil definition file", type=str)
    parser.add_argument("sim_result", help="Simulation .msh output", type=str)
    parser.add_argument("sim_geo", help="Simulation .geo output", type=str)
    parser.add_argument("--ncpus",
                        help="Number of CPUS to use in optimization",
                        type=int,
                        default=4)

    subparsers = parser.add_subparsers(help="Mode of visualization")

    direction_parser = subparsers.add_parser(
        "direction", help="Use E-Field directional optimization")
    direction_parser.set_defaults(direction=True)

    mag_parser = subparsers.add_parser(
        "magnitude", help="Use E-Field magnitude optimization")
    mag_parser.set_defaults(direction=False)

    args = parser.parse_args()
    if args.direction:
        target_direction = get_target_direction(args.json)
    else:
        target_direction = None

    msn = np.load(args.msn)

    # Run visualization in FieldOpt, simnibs is difficult..
    model = HeadModel(args.msh)
    func = FieldFunc(model,
                     None,
                     args.coil,
                     tet_weights=None,
                     direction=target_direction)
    func.visualize_evaluate(matsimnibs=msn,
                            out_sim=args.sim_result,
                            out_geo=args.sim_geo)


if __name__ == '__main__':
    main()
