#!/usr/bin/env python


from pyself.config import SelfModelConfig


## Add ctypes interface to fortran library
from ctypes import CDLL, c_int, c_double, c_char_p, POINTER, c_void_p
from ctypes.util import find_library
import os


class SelfModel:
    def __init__(self, config: SelfModelConfig, lib: str = None):
        self.case_directory = case_directory
        self.config = config

        if lib is None:
            try:
                self.lib = CDLL(find_lib("self_interface"))
            except:
                raise Exception(
                    "Could not find the libself_interface.so library. Ensure your LD_LIBRARY_PATH includes the path for libself_interface.so"
                )
        else:
            # To do - check if path exists..
            # To do - must be libself_interface.so
            self.lib = CDLL(lib)
            # To do - exceptions on  CDLL failures.

        # To do : Call fortran library to ininitialize the model
        # To do : Add mpi support - pass communicator and rank to fortran library
        InitializeModel(self.config.config_file)
        self._initialized = True

    # def report_parameters(self):

    def set_parameter(self, section: str, key: str, value: Any):
        self.config.set_parameter(section, key, value)

    def get_parameter(self, section: str, key: str, value: Any):
        return self.config.get_parameter(section, key)

    # def forward_step(self):

    # def write_pickup_file(self):

    # def get_solution(self):

    # def finalize(self):

    # def run(self):
    #     print("Running SELF model..."
