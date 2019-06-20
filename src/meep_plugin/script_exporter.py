from ipkiss3.io.vfab_export import common
from ipkiss3.io.vfab_export.common.primitives import (
            Prism,
            ObjectGroup,
            Geometry,
            Port,
            Component,
            SParamSweep,
            Difference,
            SimulationJob,
            SimulationSettings,
            Material,
            MaterialProperty
        )

import toolz as tlz
import numpy as np
import os
import ipkiss3.all as i3

import abc

class MeepExporter(object):
    __metaclass__ = abc.ABCMeta
    primitive_type = None

    def __init__(self, context):
        self.context = context

    @abc.abstractmethod
    def export(self, obj):
        pass

    def __call__(self, obj):
        return self.export(obj)

linesep = '\n'

def comment(cmt):
    return '# {}'.format(cmt)

def normalize_variable_name(s):
    """ Replaces spaces with underscores in order to
    have a 'normal' variable name.
    """
    return str.replace(s, ' ', '_')

def meep_vertices(points, z):
    """ Converts a list of 2D points into a meep vertices array.
    """
    points_str = [
        'mp.Vector3({}, {}, {})'.format(x, y, z)
        for x, y in points
    ]
    return "[" + ",".join(points_str) + "]"


class MaterialExporter(MeepExporter):
    primitive_type = Material

    def export(self, material):
        """ code to export a group of items to lumerical script"""
        context = self.context
        material_type = "Dielectric"
        material_var = normalize_variable_name(material.name)
        meep_property_map = {
                MaterialProperty.refractive_index: 'index'

        }
        self.context.add_material(material.name, material_var)
        #r, g, b, alpha = material.color;
        return [
            comment("Material '{}'".format(material.name)),
            '{} = mp.Medium({})'.format(material_var, ",".join(['{}={}'.format(meep_property_map[prop], prop_val) for prop, prop_val in material.properties.items() if prop in meep_property_map]))
        ]



class GeometryExporter(MeepExporter):
    primitive_type = Geometry

    def export(self, geometry):
        context = self.context
        for material in geometry.materials:
            for stmt in context.export(material):
                yield stmt
            yield "\n"
        _components = []
        for component in geometry.components:
            for stmt in context.export(component):
                yield stmt
            _components.append(normalize_variable_name(component.name))
            yield "\n"
        yield "_geometry=[{}]".format(",".join(_components))
        yield "geometry = [obj for g in _geometry for obj in g]"


class PrismExporter(MeepExporter):
    primitive_type = Prism
    def export(self, prism):
        """ code to export a 'prism' to meep script"""
        ctx = self.context
        z_min, z_max = prism.z_coords
        holes = list(prism.inner_points)

        # meep doesn't support defining polygons with holes directly,
        # though this can be resolved using boolean operations.

        # holes = [
        #     "{} = {};".format(hole_name(prism.name, i), lumerical_vertices(pts))
        #     for i, pts in enumerate(holes)
        # ]
        #
        # # subtract the holes from the 'exterior' polygon.
        # bool_ops = [
        #     "vertices = polydiff(vertices, {});".format(hole_name(prism.name, i))
        #     for i, _ in enumerate(holes)
        # ]

        # the temporary hole polygons can be removed.
        # remove_holes = []

        return ([
            comment('Prism {}'.format(prism.name)),
            "vertices = {}".format(meep_vertices(prism.points, z_min)),
            "{} = mp.Prism(vertices, height={}, material={}, axis=mp.Vector3(0.0, 0.0, 1.0))".format(normalize_variable_name(prism.name),
                                                                                                      z_max-z_min,
                                                                                                      ctx.get_material(prism.material))
            ]
        )


class ComponentExporter(MeepExporter):
    primitive_type = Component

    def export(self,  component):
        ctx = self.context
        yield comment("Start component {}".format(component.name))

        _objects = []
        for obj in component.objects:
            for stmt in ctx.export(obj):
                yield stmt
            yield "\n"
            _objects.append(normalize_variable_name(obj.name))
        yield "{}=[{}]".format(normalize_variable_name(component.name), ", ".join(_objects))
        yield comment("End component {}".format(component.name))


def normal_to_angles(normal):
    # returns azimuth (counterclockwise from x-axis) and inclination (from xy-plane)
    x, y, z = normal
    if np.abs(x)<1e-12:
        x = 0
    if np.abs(y)<1e-12:
        y = 0
    if np.abs(z)<1e-12:
        z = 0

    inclination = np.arcsin(z)
    if (y, z) == (0, 0):
        return np.arccos(x), 0.0
    if (x, z) == (0, 0):
        return np.arcsin(y), 0.0
    if (x, y) == (0, 0):
        return 0.0, inclination
    if x == 0.0:
        return 0.5 * np.pi * np.sign(y), inclination
    return np.arctan2(y , x), inclination


class PortExporter(MeepExporter):
    primitive_type = Port

    def export(self, port):
        normal = [
            0.0 if np.isclose(n, 0.0)
            else n
            for n in port.normal
        ]
        if np.count_nonzero(normal) > 1:
            raise ValueError("Port `{}`'s orientation is not supported, please change its normal to be parallel to the calculation domain's boundaries, choose between X, Y or Z.".format(port.name))

        ctx = self.context
        x, y, z = port.position
        azimuth, inclination = normal_to_angles(port.normal)

        port_width, port_height = port.box_size
        port_name = port.name

        # only ports with inclination 0 and 90 supported so far, or rather ports with normal along either X,Y or Z
        if np.isclose(inclination, 0.0):
            mwx = abs(port_width * np.sin(azimuth))
            mwy = abs(port_width * np.cos(azimuth))
            mwz = abs(port_height)
        elif np.isclose(abs(90.0-inclination), 0.0):
            mwz = 0.0
            mwx = abs(port_height * np.cos(azimuth)) + abs(port_width * np.sin(azimuth))
            mwy = abs(port_width * np.cos(azimuth)) + abs(port_height * np.sin(azimuth))
        else:
            raise ValueError("Free ports not supported yet")

        if np.isclose(mwx, 0.0):
            mwx = 0
            if np.cos(azimuth) < 0:
                orientation = 'xmin'
            else:
                orientation = 'xmax'
        if np.isclose(mwy, 0.0):
            mwy = 0.0
            if np.sin(azimuth) < 0:
                orientation = 'ymin'
            else:
                orientation = 'ymax'
        if np.isclose(mwz, 0.0):
            mwz = 0.0
            if np.sin(inclination) < 0:
                orientation = 'zmin'
            else:
                orientation = 'zmax'

        xmin, xmax = [x - mwx / 2.0, x + mwx / 2.0]
        ymin, ymax = [y - mwy / 2.0, y + mwy / 2.0]
        zmin, zmax = [z - mwz / 2.0, z + mwz / 2.0]

        wlcen = 1.55
        fcenter = 1 / wlcen 

        params = {
            'port_name': port_name,
            'mwx': mwx,
            'x': x,
            'y': y, 
            'z': z,
            'mwz': mwz,
            'mwy': mwy,
            'fcenter': fcenter,
            'fwidth': 0.2 * fcenter
        }

        yield "port_size_{port_name} = mp.Vector3({mwx}, {mwy}, {mwz})".format(**params)
        yield "port_center_{port_name} = mp.Vector3({x}, {y}, {z})".format(**params)
        #yield "port_volume_{port_name} = mp.Block(port_size_{port_name}, center=port_center_{port_name})".format(**params)
#        yield "port_monitor_{port_name} = simulation.add_mode_monitor({fcenter}, 0, 1, mp.ModeRegion(volume=port_volume_{port_name}))".format(**params)

        yield """
port_src_{port_name} = mp.EigenModeSource(src=mp.GaussianSource({fcenter}, fwidth={fwidth}),
    size=port_size_{port_name},
    center=port_center_{port_name},
    eig_band=1,
    eig_parity=mp.NO_PARITY,
    eig_match_freq=True
)""".format(**params)
 

class SParamSweepExporter(MeepExporter):
    primitive_type = SParamSweep

    def export(self, ssweep):
        ctx = self.context

        yield "# TODO: implement sparam sweep" 
        yield "#ports_to_sweep = port_monitors.keys()"
        yield "# from meep_plugin.utils import sparam_sweep"
        yield "# sweep_result = sparam_sweep(simulation, port_monitors, port_sources, ports_to_sweep)" 



class SimulationSettingsExporter(MeepExporter):
    primitive_type = SimulationSettings
    def export(self, simsettings):
        ctx = self.context
        #
        start_wl = 1.5e-6
        # stop_wl = 1.6e-6
        # mesh_accuracy = 2
        (xbox_min, xbox_max), (ybox_min, ybox_max), (zbox_min, zbox_max) = simsettings.bounding_box
        yield "bbox = ((xbox_min, xbox_max), (ybox_min, ybox_max), (zbox_min, zbox_max)) = (({}, {}), ({}, {}), ({}, {}))".format(xbox_min, xbox_max, ybox_min, ybox_max, zbox_min, zbox_max)
        resolution = 40
        yield "cell = mp.Vector3({}, {}, {})".format(xbox_max - xbox_min, ybox_max - ybox_min, zbox_max - zbox_min)
        yield """
geometry_center = mp.Vector3(
    xbox_min + (xbox_max - xbox_min) / 2,
    ybox_min + (ybox_max - ybox_min) / 2,
    zbox_min + (zbox_max - zbox_min) / 2
)
"""

        sim_setup = """simulation = mp.Simulation(cell_size=cell,
        boundary_layers=[mp.PML(.2)],
        sources=[port_src_{port_src}],
        resolution={resolution},
        geometry=geometry,
        geometry_center=geometry_center
        )
        """.format(port_src=simsettings.ports[0].name, resolution=resolution)

        yield sim_setup
        # sim_setup = simulation_setup(xbox, ybox, zbox,
        #                              start_wl, stop_wl,
        #                              mesh_accuracy=mesh_accuracy)
        yield "simulation.init_sim()"

        yield """

def _mk_pos_convertor(sim, bbox):
    eps_data = sim.get_array(component=mp.Dielectric)
    shp = eps_data.shape
    (xbox_min, xbox_max), (ybox_min, ybox_max), (zbox_min, zbox_max) = bbox

    x_cell = (xbox_max - xbox_min) / shp[0]
    y_cell = (ybox_max - ybox_min) / shp[1]
    z_cell = (zbox_max - zbox_min) / shp[2]

    def to_cell_pos(x, y, z):
        import numpy as np
        if x < xbox_min or x > xbox_max:
            x_idx = None
        else:
            x_rel = x - xbox_min
            x_idx = np.round(x_rel / x_cell)

        if y < ybox_min or y > ybox_max:
            y_idx = None
        else:
            y_rel = y - ybox_min
            y_idx = np.round(y_rel / y_cell)
 
        if z < zbox_min or z > zbox_max:
            z_idx = None
        else:
            z_rel = z - zbox_min
            z_idx = np.round(z_rel / z_cell)

        return (int(x_idx), int(y_idx), int(z_idx))

    return to_cell_pos 

to_cell_pos = _mk_pos_convertor(simulation, bbox)
"""



class SimulationJobExporter(MeepExporter):
    primitive_type = SimulationJob

    def export(self, simjob):
        ctx = self.context

        yield 'import meep as mp'
        yield 'import matplotlib.pyplot as plt'
        yield 'import numpy as np'

        for stmt in ctx.export(simjob.geometry):
            yield stmt

        yield '\n'


        simsettings = simjob.settings
        for port in simsettings.ports:
            for stmt in ctx.export(port):
                yield stmt

        port_names = [port.name for port in simsettings.ports]

        yield "port_names = {}".format(port_names)

        yield """
port_sources = {
    port_name: locals().get('port_src_' + port_name)
    for port_name in port_names
}
"""
        yield """
port_volumes = {
    port_name: locals().get('port_volume_' + port_name)
    for port_name in port_names
}
"""
        for stmt in ctx.export(simjob.settings):
            yield stmt

        for sweep in simsettings.sweeps:
            for stmt in ctx.export(sweep):
                yield stmt

        for macro in simsettings.macros:
            for stmt in macro.commands:
                yield stmt


def get_export_map(context):
    ExportClasses = MeepExporter.__subclasses__()
    return {
        exportcls.primitive_type: exportcls(context)
        for exportcls in ExportClasses
    }


class Context(object):

    def __init__(self, top, export_map=None):
        self.export_map = export_map or get_export_map(self)
        self.top = top
        self._materials = {}

    def export(self, inst):
        inst_type = type(inst)

        export_fn = self.export_map[inst_type]
        return export_fn(inst)

    def get_identifier(self, obj):
        pass

    def get_material(self, name):
        return self._materials.get(name, None)

    def add_material(self, name, variable):
        self._materials[name] = variable

def export(geometry):
    ctx = Context(top=geometry)
    lines = list(ctx.export(geometry))
    return lines
