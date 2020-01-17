import numpy as np

import clipper_python

from mdc3.types.real_space import MCDXMap


def ccp4_path_to_np(ccp4_path):
    xmap = load_ccp4_map(ccp4_path)

    xmap_np = xmap.export_numpy()

    return xmap_np


def make_mean_map(xmaps_np):
    # values = np.array([xmap.get_map_data(sparse=False).as_numpy_array() for xmap in xmaps])
    # print(values.shape)
    # mean_values = np.mean(values, axis=0)
    # mean_map = xmaps[0].new_from_template(map_data=flex.double(mean_values))

    mean_map_np = np.mean(xmaps_np, axis=0)

    return mean_map_np


def output_mean_map(template_xmap, mean_map_np, path):
    print("Outputting mean map to: {}".format(path))
    spacegroup = template_xmap.xmap.spacegroup
    cell = template_xmap.xmap.cell
    grid = template_xmap.xmap.grid_sampling
    grid_params = (grid.nu,
                   grid.nv,
                   grid.nw,
                   )
    # resolution = template_xmap.hkl_info.resolution
    mean_map_xmap = MCDXMap.xmap_from_numpy(spacegroup,
                                            cell,
                                            # resolution,
                                            grid_params,
                                            mean_map_np,
                                            )

    mean_map_xmap.to_ccp4(path)


def output_mean_nxmap(mean_map_np, cell, path, grid_params):
    print("Outputting mean map to: {}".format(path))
    rot = clipper_python.Mat33_double(np.eye(3))
    trans = clipper_python.Vec3_double(0.0, 0.0, 0.0)

    rtop = clipper_python.RTop_orth(rot,
                                    trans,
                                    )

    # Generate the clipper grid
    grid = clipper_python.Grid(grid_params[0],
                               grid_params[1],
                               grid_params[2],
                               )
    # resolution = template_xmap.hkl_info.resolution
    mean_map_nxmap = clipper_python.NXmap_float(grid,
                                                rtop,
                                                )

    mean_map_nxmap.import_numpy(clipper_python.Coord_grid(0, 0, 0),
                                mean_map_np,
                                )

    mapout = clipper_python.CCP4MAPfile()
    # print(dir(mapout))
    mapout.open_write(str(path))
    mapout.set_cell(cell)
    mapout.export_nxmap_float(mean_map_nxmap)
    mapout.close_write()


def load_ccp4_map(ccp4_path):
    xmap = clipper_python.Xmap_float()

    mapout = clipper_python.CCP4MAPfile()
    mapout.open_read(str(ccp4_path))
    mapout.import_xmap_float(xmap)
    mapout.close_read()
    return xmap


def save_nxmap_from_template(template_xmap, mean_map_np, path):
    print("Outputting mean map to: {}".format(path))
    cell = template_xmap.cell
    grid = template_xmap.grid_sampling
    grid_params = [grid.nu,
                   grid.nv,
                   grid.nw,
                   ]
    rot = clipper_python.Mat33_double(np.eye(3))
    trans = clipper_python.Vec3_double(0.0, 0.0, 0.0)

    rtop = clipper_python.RTop_orth(rot,
                                    trans,
                                    )

    # Generate the clipper grid
    grid = clipper_python.Grid(grid_params[0],
                               grid_params[1],
                               grid_params[2],
                               )
    # resolution = template_xmap.hkl_info.resolution
    mean_map_nxmap = clipper_python.NXmap_float(grid,
                                                rtop,
                                                )

    mean_map_nxmap.import_numpy(clipper_python.Coord_grid(0, 0, 0),
                                mean_map_np,
                                )

    mapout = clipper_python.CCP4MAPfile()
    # print(dir(mapout))
    mapout.open_write(str(path))
    mapout.set_cell(cell)
    mapout.export_nxmap_float(mean_map_nxmap)
    mapout.close_write()