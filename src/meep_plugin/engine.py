
from ipkiss3.simulation.device import EMSimulation
from ipkiss3.simulation.device._engine import EMSimulationEngine



class MEEPEngine(EMSimulationEngine):

    def inspect(self, simdata):
        """ Inspects the simulation job in the tool.
        """
        from .script_exporter import export
        simulation = export(simdata)
        return simulation

    def get_result(self, simdata):
        raise NotImplementedError()
        

class MEEPFDTDSimulation(EMSimulation):

    def _get_simdata(self):
        from ipkiss3.simulation.device._export import _export_simulation
        sim_data = _export_simulation(self, holes_support=False)
        return sim_data

    def _default_engine(self):
        return MEEPEngine()

    def inspect(self):
        """ Inspect the simulation job in the tool.
        """
        simdata = self._get_simdata()
        return self.engine.inspect(simdata)



