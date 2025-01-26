from addon import OpenmmAddOn
from gi.repository import AmlCore, AmlLammpsIo
from pathlib import Path

class DumpActionParams(AmlCore.ActionParams):
    def __init__(self, binary=True):
        super().__init__()
        self.action = AmlLammpsIo.BinaryDumpWriter() if binary else AmlLammpsIo.DumpWriter()
        self.action_params = AmlLammpsIo.WriterParams()
        self.dumpdir = ""

    def do_copy(self):
        res = DumpActionParams()
        res.action = self.action
        res.action_params = self.action_params
        res.dumpdir = self.dumpdir
        return res

class DumpAction(AmlCore.Action):
    def __init__(self):
        super().__init__()

    def do_get_params_error_message(self, params):
        return ""

    def do_perform(self, data : AmlCore.DataCollection):
        params = self.get_params()
        writer = params.action
        wp = params.action_params

        dumpdir_path = Path(data.get_element("filesystem.root").get_val()).joinpath(params.dumpdir)
        dumpdir_path.mkdir(parents=True, exist_ok=True)
        step = data.get_element("repr.step").get_val()
        wp.set_filepath(str(dumpdir_path.joinpath(f"{step}.lammpsdump")))
        
        writer.set_params(wp)

        writer.perform(data)

class LammpsDump(OpenmmAddOn):
    def __init__(self):
        self.action = DumpAction()
        self.params = None

    def init_dict(self, dictionary):
        self.params = DumpActionParams(dictionary.get("binary", True))
        self.params.dumpdir = dictionary["directory"]
        self.params.action_params.set_box_id("repr.box")
        self.params.action_params.set_timestep_id("repr.step")
        self.params.action_params.set_particles_id("repr.particles")
        self.params.action_params.set_properties(dictionary["columns"])
        self.action.set_params(self.params)

    def get_action(self):
        return self.action
