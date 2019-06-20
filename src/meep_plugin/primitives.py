from ipkiss3.io.vfab_export.common.primitives import SimulationSettings
from atom import api


class MeepSimulationSettings(SimulationSettings):
    resolution = api.Typed(float)