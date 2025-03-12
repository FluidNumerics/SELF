#!/usr/bin/env python


from pyself.config import SelfModelConfig


## Add ctypes interface to fortran library
from ctypes import CDLL, c_int, c_double, c_char_p, POINTER, c_void_p
from ctypes.util import find_library
import os


class SelfModel:
    def __init__(self, config: SelfModelConfig = SelfModelConfig(), lib: str = None):
        self.case_directory = case_directory
        self.config = config
        self._config_file = f"{self.config.case_directory}/model_input.json"

        if lib is None:
            try:
                self._lib = CDLL(find_lib("self_interface"))
            except:
                raise Exception(
                    "Could not find the libself_interface.so library. Ensure your LD_LIBRARY_PATH includes the path for libself_interface.so"
                )
        else:
            # Library ust be libself_interface.so
            if not lib.endswith("libself_interface.so"):
                raise Exception("Library must be libself_interface.so")
            # Library must exist
            if not os.path.exists(lib):
                raise Exception(f"Could not find the library {lib}")
            try:
                self._lib = CDLL(lib)
            except:
                raise Exception(f"Could not load the library {lib}")

        self._configure_interface()

        self._initialized = False

    def _configure_interface(self):
        self._lib.Initialize.argtypes = [c_char_p]
        self._lib.Initialize.restype = None

        self._lib.UpdateParameters.argtypes = []
        self._lib.UpdateParameters.restype = None

        self._lib.ForwardStep.argtypes = [c_double, c_double]  # No arguments
        self._lib.ForwardStep.restype = c_int  # Function returns an integer

        self.lib.WritePickupFile.argtypes = [c_char_p]
        self.lib.WritePickupFile.restype = None

        # self.lib.get_solution.argtypes = [c_void_p]
        # self.lib.get_solution.restype = POINTER(c_double)

        # self.lib.finalize.argtypes = [c_void_p]
        # self.lib.finalize.restype = None

    # def report_parameters(self):

    def set_parameter(self, section: str, key: str, value: Any):
        self.config.set_parameter(section, key, value)

    def get_parameter(self, section: str, key: str, value: Any):
        return self.config.get_parameter(section, key)

    def update_parameters(self):
        self.config.save_config()
        if self._initialized:
            self._lib.UpdateParameters()
        else:
            raise Exception(
                "Configuration file saved, but not pushed to model. Model is not initialized"
            )

    def initialize_model(self):
        if not self._initialized:
            # Save the config to the case directory
            self.config.save_config()
            # Call the initialize model function
            self._lib.Initialize(self._config_file.encode("utf-8"))

            # To do, print out model parameters, nicely formatted
            self._initialized = True
        else:
            raise Exception("Model is already initialized")

    def forward_step(self, dt, update_interval):
        if not self._initialized:
            self.initialize_model()

        err = self._lib.ForwardStep(c_double(dt), c_double(update_interval))

        # To do:  error handling

    # def write_pickup_file(self):

    # def get_solution(self):

    # def finalize(self):

    # def run(self):
    #     print("Running SELF model..."
