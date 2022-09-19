from collections import namedtuple
from enum import Enum

import argparse
import json
import numpy as np

import gmsh
import vtk

import pyvista as pv
pv.start_xvfb()

ENORM = 1
WEIGHTFUNC = 2

Pole = namedtuple("Pole", ["start", "end"])
CoilStick = namedtuple("CoilStick", ["arrow", "pole"])
CameraPos = namedtuple("CameraPos", ["base", "direction"])


def get_field_data(elm_ids):

    _, norm_inds, norms, *_ = gmsh.view.get_homogeneous_model_data(tag=ENORM,
                                                                   step=0)
    _, wf_inds, wf, *_ = gmsh.view.get_homogeneous_model_data(tag=WEIGHTFUNC,
                                                              step=0)

    return (norms[elm_ids], wf[elm_ids])


def create_brain_model(mesh_file):

    gmsh.initialize()
    gmsh.open(mesh_file)

    node_ids, node_coords, node_params = gmsh.model.mesh.get_nodes(tag=2)
    node_coords = node_coords.reshape(-1, 3)
    _, elm_ids, elm_node_ids = gmsh.model.mesh.get_elements(tag=2, dim=3)
    elm_ids = elm_ids[0]
    elm_node_ids = elm_node_ids[0]

    points = node_coords

    n_cells = elm_ids.shape[0]
    cells = np.empty((n_cells, 5), dtype=int)

    cells[:, 0] = 4
    cells[:, 1:] = elm_node_ids.reshape(-1, 4) - 1
    cells = cells.ravel()

    celltypes = np.empty(n_cells, dtype=np.uint8)
    celltypes[:] = vtk.VTK_TETRA

    grid = pv.UnstructuredGrid(cells, celltypes, points)
    norms, wf = get_field_data([elm_ids - 1])

    gmsh.clear()
    gmsh.finalize()
    return grid, norms, wf


def make_coil_stick(orientation):

    centre = orientation[:3, -1]
    coil = orientation[:3, 1]
    normal = orientation[:3, 2]

    return CoilStick(arrow=pv.Arrow(start=centre, direction=coil, scale=10),
                     pole=pv.Tube(pointa=centre,
                                  pointb=centre + (normal * 50),
                                  radius=0.25))


def get_directional_arrow(coordinate, direction):
    return pv.Arrow(start=coordinate, direction=direction, scale=10)


def add_stick_to_plotter(plotter, stick, **kwargs):

    plotter.add_mesh(stick.arrow, **kwargs)
    plotter.add_mesh(stick.pole, **kwargs)
    return


def camera_from_msn(orientation):
    base = orientation[:3, -1]
    direction = -orientation[:3, 2]

    return CameraPos(base, direction)


class OptType(Enum):
    magnitude = "magnitude"
    direction = "direction"

    def __str__(self):
        return self.value


def position_and_direction_from_spec(json_file):

    with open(json_file, 'r') as f:
        specs = json.load(f)

    coordinate = np.array([specs['pos_x'], specs['pos_y'], specs['pos_z']],
                          dtype=float)

    direction = np.array([specs['dir_x'], specs['dir_y'], specs['dir_z']],
                         dtype=float)

    return coordinate, direction


def main():
    p = argparse.ArgumentParser(description="Generate a QC visualization"
                                " of the ADM optimization method")
    p.add_argument("sim_msh", help="Simulation mesh file", type=str)
    p.add_argument("matsimnibs", help="Optimal matsimnibs matrix", type=str)
    p.add_argument("sim_spec_json",
                   help="Path to simulation specification JSON file",
                   type=str)
    p.add_argument("--export-img", help="QC image to output")
    p.add_argument("--export-html", help="Interactive HTML to output")
    p.add_argument("optimization_type", type=OptType, choices=list(OptType))

    args = p.parse_args()

    msn = np.load(args.matsimnibs)
    coord, direction = position_and_direction_from_spec(args.sim_spec_json)

    coilstick = make_coil_stick(msn)
    camera = camera_from_msn(msn)

    grid, norms, wf = create_brain_model(args.sim_msh)
    grid.cell_data['norms (V/m)'] = norms
    # grid.cell_data['target'] = wf

    p = pv.Plotter(polygon_smoothing=True, off_screen=True)
    p.add_mesh(grid, cmap="jet", smooth_shading=True)

    add_stick_to_plotter(p, coilstick, color="green")

    additional_text = ""
    if args.optimization_type == OptType.direction:
        direction_arrow = get_directional_arrow(msn[:3, -1], direction)
        normal_arrow = get_directional_arrow(coord, -direction)
        p.add_mesh(direction_arrow, color="black")
        p.add_mesh(normal_arrow, color="white")
        additional_text = (
            "\nBlack arrow indicates field target direction"
            "\nWhite arrow shows computed target with computed normal vector"
            "\nMagnitude shown is magnitude along"
            " the direction of the black arrow")

    p.camera.position = camera.base + (camera.direction * 100)
    p.enable_anti_aliasing()
    p.window_size = [1000, 1000]
    p.camera.zoom(0.50)
    p.add_text("Green arrow indicates coil anterior position"
               "\nScalar map is E-field magnitude (V/m)" + additional_text,
               color="white",
               shadow=True,
               font_size=15)

    if args.export_img:
        p.screenshot(args.export_img)

    if args.export_html:
        p.export_html(args.export_html)


if __name__ == '__main__':
    main()
