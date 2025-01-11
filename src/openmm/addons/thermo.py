from addon import OpenmmAddOn, OpenmmSimulation
from gi.repository import AmlCore
from pathlib import Path
import os

class ThermoActionParams(AmlCore.ActionParams):
    def __init__(self):
        super().__init__()

    def do_copy(self):
        return ThermoActionParams()

class ThermoAction(AmlCore.Action):
    def __init__(self):
        super().__init__()

    def do_get_params_error_message(self, params):
        return ""

    def do_perform(self, data : AmlCore.DataCollection):
        step = data.get_element("repr.step").get_val()
        temperature = data.get_element("repr.thermo.temperature").get_val()
        kinetic = data.get_element("repr.thermo.kinetic").get_val()
        potential = data.get_element("repr.thermo.potential").get_val()
        total = data.get_element("repr.thermo.total").get_val()

        thermofile = data.get_element("filesystem.thermofile").get_val()

        if not os.path.isfile(thermofile):
            with open(thermofile, "w") as f:
                f.write("step,potential,kinetic,total,temperature\n")

        with open(thermofile, "a") as f:
            f.write(f"{step},{potential},{kinetic},{total},{temperature}\n")


class ThermoDump(OpenmmAddOn):
    def __init__(self):
        self.action = ThermoAction()
        self.params = ThermoActionParams()

    def init_dict(self, dictionary):
        self.action.set_params(self.params)

    def get_action(self):
        return self.action
