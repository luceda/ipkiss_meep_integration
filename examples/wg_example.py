import ipkiss3.all as i3
import numpy as np
import os


wg_process = i3.ProcessLayer(name="wg", extension="WG")
core_purpose = i3.PatternPurpose(name="core", extension="CORE")
clad_purpose = i3.PatternPurpose(name="core", extension="CLAD")
core_layer = i3.ProcessPurposeLayer(process=wg_process,
                                    purpose=core_purpose)
clad_layer = i3.ProcessPurposeLayer(process=wg_process,
                                    purpose=clad_purpose)




def vfab():
    from ipkiss.plugins.vfabrication.process_flow import VFabricationProcessFlow
    from pysics.basics.material.material import Material, MaterialFactory
    from pysics.basics.material.material_stack import MaterialStackFactory, MaterialStack
    from ipkiss.visualisation.display_style import DisplayStyle
    from ipkiss.visualisation.color import COLOR_BLUE, COLOR_RED
    wg_layer = clad_layer - core_layer

    silicon = Material(name="silicon", epsilon=3.5 ** 2, display_style=DisplayStyle(color=COLOR_BLUE))
    oxide = Material(name="silicon oxide", epsilon=1.5 ** 2, display_style=DisplayStyle(color=COLOR_RED))

    silicon = MaterialStack(name="silicon",
                            materials_heights=[
                                (oxide, 0.5),
                                (silicon, 0.2),
                                (oxide, 0.5)
                            ])

    oxide = MaterialStack(name="oxide",
                          materials_heights=[
                              (oxide, 0.5),
                              (oxide, 0.2),
                              (oxide, 0.5)
                          ])

    vfab_process = VFabricationProcessFlow(active_processes=[wg_process],
                                           process_layer_map={wg_process: wg_layer},
                                           process_to_material_stack_map= # (wg_process, )
                                           [((0,), silicon),
                                            ((1,), oxide),
                                           ],
                                           is_lf_fabrication={wg_process: False}
                                          )
    return vfab_process



wg_width = 0.5
wg_length = 2.0

# waveguide template
wg_tmpl = i3.WindowWaveguideTemplate()
wg_tpl_lay = wg_tmpl.Layout(windows=[
    i3.PathTraceWindow(layer=core_layer,
                       start_offset=-0.5 * wg_width,
                       end_offset=+0.5 * wg_width),
    i3.PathTraceWindow(layer=clad_layer,
                       start_offset=-0.5 * wg_width - 2.0,
                       end_offset=0.5 * wg_width + 2.0)
    ]
)

vfab_process = vfab()
xsection = wg_tpl_lay.cross_section(process_flow=vfab_process)
xsection.visualize()

# waveguide
wg = i3.Waveguide(trace_template=wg_tmpl)
wg_lo = wg.Layout(shape=[(-wg_length /2., 0.0), (wg_length/2., 0.0)])
wg_lo.visualize_2d(vfabrication_process_flow=vfab_process)


wg_lo_si = wg_lo.size_info()
geometry = i3.device_sim.SimulationGeometry(
    layout=wg_lo,
    process_flow=vfab_process,
    bounding_box=[
        (wg_lo_si.west - 1.0, wg_lo_si.east + 1.0),
        (wg_lo_si.south - 1.0, wg_lo_si.north + 1.0),
        None,
    ]
)

from meep_plugin import MEEPFDTDSimulation

simjob = MEEPFDTDSimulation(
    geometry=geometry,
    monitors=[i3.device_sim.Port(name="in", n_modes=2),
              i3.device_sim.Port(name="out", n_modes=2)],
    outputs=[
        i3.device_sim.SMatrixOutput(
            name='smatrix',
            wavelength_range=(1.50, 1.55, 20),
            symmetries=[('out', 'in')]
        )
        ]
)

stmts = simjob.inspect()
#sweep = ctx.sweeps['sparam']

with open('meep_script.py', 'w') as out:
    out.write('\n'.join(stmts))

"""
smat = simjob.get_result(name='smatrix')
# Plot result
wavelengths = smat.sweep_parameter_values
dB = lambda x: 20*np.log10(np.abs(x))

import matplotlib.pyplot as plt
plt.figure()
plt.plot(wavelengths, dB(smat["in", "out"]), label='transmission', linewidth=5)
plt.plot(wavelengths, dB(smat["in", "in"]), label='reflection', linewidth=5)
plt.legend()
plt.xlabel('wavelength [$\mu m$]')
plt.ylabel('S-parameter magnitude')
plt.show()
"""

