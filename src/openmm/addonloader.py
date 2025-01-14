from pathlib import Path
import importlib
import importlib.util
import sys
import inspect
from addon import OpenmmAddOn


class AddOnLoader:
    def __init__(self, logger=None):
        self.logger = logger

        if self.logger:
            self.logger.info("Initializing AddOnLoader")

        self.dir = Path(__file__).parent.joinpath("addons")
        self.addon_classes = {}

    def load_addons_from_file(self, path):
        if self.logger:
            self.logger.info(f"Processing file \"{path}\"")

        spec = importlib.util.spec_from_file_location(path.name, path)
        foo = importlib.util.module_from_spec(spec)
        sys.modules[path.name] = foo
        spec.loader.exec_module(foo)

        for name, obj in inspect.getmembers(foo):
            if inspect.isclass(obj) and issubclass(obj, OpenmmAddOn) and name not in self.addon_classes:
                if self.logger:
                    self.logger.info(f"Found AddOn class \"{name}\"")
                self.addon_classes[name] = obj

    def load(self):
        sources = self.dir.glob("*.py")
        
        if self.logger:
            self.logger.info(f"Loading addons. Found {len(list(sources))} files")

        for path in sources:
            self.load_addons_from_file(path)

    def get_addon_class(self, name):
        return self.addon_classes[name]
