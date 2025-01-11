# from addon import OpenmmAddOn, OpenmmSimulation
# from gi.repository import AmlCore
# from pathlib import Path

# class CheckpointActionParams(AmlCore.ActionParams):
#     def __init__(self):
#         pass

# class CheckpointAction(AmlCore.Action):
#     def get_params_error_message(self, params):
#         return ""

#     def perform(self, data : AmlCore.DataCollection):
#         openmm_object = data.get_element("openmm")
#         checkpointdir_path = Path(data.get_element("filesystem.checkpointdir").get_val())
#         step = data.get_element("repr.step").get_val()

#         with open(checkpointdir_path.joinpath(f"{step}.openmmcheckpoint"), "wb") as f:
#             f.write(openmm_object.context.createCheckpoint())

# class MakeCheckpoint(OpenmmAddOn):
#     def __init__(self):
#         self.action = CheckpointAction()
#         self.params = CheckpointActionParams()

#     def init_dict(self, dictionary):
#         self.action.set_params(self.params)

#     def get_action(self):
#         return self.action
