#!/usr/bin/env python


from pyself.config import SelfModelConfig


## Add ctypes interface to fortran library
from ctypes import CDLL, c_int, c_double, c_char_p, POINTER, c_void_p
import os

# Load the fortran library
lib = CDLL(os.path.join(os.path.dirname(__file__), "libself_interface.so"))


class SelfModel:
    def __init__(self, config: SelfModelConfig):
        self.case_directory = case_directory
        self.config = config

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
    #     print("Running SELF model...")
