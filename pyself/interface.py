#!/usr/bin/env python


from pyself.config import SelfModelConfig


## Add ctypes interface to fortran library
from ctypes import (
    CDLL,
    c_int,
    c_double,
    c_char_p,
    POINTER,
    c_void_p,
    create_string_buffer,
)
import numpy as np
from ctypes.util import find_library
import os

_VAR_BUFFER_SIZE = 256


class SelfModel:
    def __init__(self, config: SelfModelConfig = SelfModelConfig(), lib: str = None):
        self.case_directory = case_directory
        self.config = config
        self._config_file = f"{self.config.case_directory}/model_input.json"
        self._solution = None
        self._mesh = None
        self._geometry = None
        self._last_pickup_file = None

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
        self._precision = self._lib.GetPrecision()
        # self._dtype = {4: np.float32, 8: np.float64}[self._precision]
        # self._cprec = {4: c_float, 8: c_double}[self._precision]

        self._initialized = False

    def _configure_interface(self):
        """Private method to configure the interface to the Fortran library"""

        self._lib.Initialize.argtypes = [c_char_p]
        self._lib.Initialize.restype = c_int

        self._lib.UpdateParameters.argtypes = []
        self._lib.UpdateParameters.restype = None

        self._lib.ForwardStep.argtypes = [c_double, c_double]  # No arguments
        self._lib.ForwardStep.restype = c_int  # Function returns an integer

        self._lib.WritePickupFile.argtypes = [c_char_p, c_char_p]
        self._lib.WritePickupFile.restype = None

        self._lib.GetSolution.argtypes = [
            POINTER(c_void_p),
            POINTER(c_int * 5),
            POINTER(c_int),
        ]
        self._lib.GetSolution.restype = None  # Subroutine, no return

        self._lib.GetVariableName.argtypes = [c_int, c_char_p]
        self._lib.GetVariableName.restype = None

        self._lib.GetPrecision.argtypes = []  # No arguments
        self._lib.GetPrecision.restype = c_int  # Function returns an integer

        self.lib.Finalize.argtypes = []
        self.lib.Finalize.restype = None

    def report_config(self):
        """Print the configuration to the console."""
        print("=" * 40)
        print(" Model Configuration ".center(40, "="))
        print("=" * 40)

        print("\n[Model]")
        model_name = self.config.config["model_name"]
        print(f"  SELF Configuration Version : {self.config.config['version']}")
        print(f"  Model Name : {model_name}")
        print(f"  Case Directory : {self.config.case_directory}")
        print(f"  Config File : {self._config_file}")
        print(f"  Precision : {self._precision}")
        print(f"  Initialized : {self._initialized}")

        print("\n[Geometry]")
        for key, value in self.config.config["time_options"].items():
            print(f"  {key.replace('_', ' ').capitalize()} : {value}")

        print("\n[Time Options]")
        for key, value in self.config.config["time_options"].items():
            print(f"  {key.replace('_', ' ').capitalize()} : {value}")

        print("\n[{model_name}]")
        for key, value in self.config.config[model_name].items():
            print(json.dumps(value, indent=4))

    def set_parameter(self, section: str, key: str, value: Any):
        """Set a specific parameter within the model configuration."""
        self.config.set_parameter(section, key, value)

    def get_parameter(self, section: str, key: str, value: Any):
        """Retrieve a specific parameter from the model configuration."""
        return self.config.get_parameter(section, key)

    def update_parameters(self):
        """Push the configuration to the model by writing the configuration
        to file and calling the Fortran-side UpdateParameters function."""

        self.config.save_config()
        if self._initialized:
            self._lib.UpdateParameters()
        else:
            raise Exception(
                "Configuration file saved, but not pushed to model. Model is not initialized"
            )

    def set_time_integrator(self, integrator: str):
        """Set the time integrator in the configuration file. The
        selfModel.config attribut is updated and saved to the case directory
        json file.

        Parameters:
        -----------
        integrator (str): the time integrator to use. Must be one of
                          'euler', 'rk2', 'rk3', or 'rk4'
        """

        self.config.set_parameter("time_options", "integrator", integrator)
        self.config.save_config()

    def initialize_model(self):
        """Initialize the model by calling the Fortran Initialize function.
        On the Fortran side this allocates memory and sets up the appropriate
        data structures for the model, based on the configuration file.
        On exit, the self._initialized flag is set to True."""

        if not self._initialized:
            # Save the config to the case directory
            self.config.save_config()
            # Call the initialize model function
            error = self._lib.Initialize(self._config_file.encode("utf-8"))
            if error != 0:
                raise Exception(
                    f"Model returned error code {error} for model_name = {self.config.config['model_name']}"
                )

            # To do, print out model parameters, nicely formatted
            self._initialized = True
        else:
            raise Exception("Model is already initialized")

    def finalize_model(self):
        """Finalize the model by calling the Fortran Finalize function.
        On the Fortran side this deallocates memory and cleans up the model.
        On exit, the self._initialized flag is set to False."""

        if self._initialized:
            self._lib.Finalize()
            self._initialized = False
        else:
            raise Exception("Model is not initialized")

    def forward_step(self, dt, update_interval):
        """Advance the model forward in time by calling the Fortran ForwardStep function.
        The function takes two arguments: the time step dt and the number of time steps to take.
        The function returns an error code, which is 0 if the function executed successfully.
        The time integrator is controlled by the configuration file in the
        time_options.integrator setting"""

        if not self._initialized:
            self.initialize_model()

        if self._precision == 4:
            err = self._lib.ForwardStep(c_float(dt), c_float(update_interval))
        else:
            err = self._lib.ForwardStep(c_double(dt), c_double(update_interval))

        return err

    def write_pickup_file(self):
        """Write the pickup file by calling the Fortran WritePickupFile function.
        The function takes a case directory as an argument and returns the name of the pickup file.
        The pickup file is written to the case directory and follows the format "solution.X.h5", where
        "X" is a 13-digit zero padded integer that corresponds to the iterate number in the simulation.

        Returns:
        --------
        pickup_file (str): the name of the pickup file that was written to disk.

        """
        if not self._initialized:
            raise Exception("Model is not initialized")

        # Prepare input and output buffers
        case_directory = self.config.case_directory.encode(
            "utf-8"
        )  # Convert string to bytes (null-terminated)
        pickup_file_buffer = create_string_buffer(buffer_size)  # Preallocated buffer

        # Call the Fortran subroutine
        self._lib.WritePickupFile(case_directory, pickup_file_buffer)

        pickup_file = pickup_file_buffer.value.decode("utf-8").strip()
        self._last_pickup_file = pickup_file

        # Convert to Python string (remove null terminator and spaces)
        return pickup_file

    def get_solution(self):
        if not self._initialized:
            raise Exception("Model is not initialized")

        # Get the solution
        solution_ptr = c_void_p()
        shape = (c_int * 5)()
        rank = c_int()
        precision = c_int()
        self._lib.GetSolution(byref(solution_ptr), shape, byref(rank))

        # Extract shape values
        dim = [shape[i] for i in range(rank.value)]  # Extract only relevant dimensions

        # Convert void pointer to float or double pointer
        if self._precision == 4:
            data_ptr = ctypes.cast(solution_ptr, POINTER(c_float))
        else:
            data_ptr = ctypes.cast(solution_ptr, POINTER(c_double))

        # Convert to NumPy array (handling column-major storage)
        solution = np.ctypeslib.as_array(
            data_ptr, shape=tuple(reversed(dim))
        )  # Reverse shape for row-major order

        self._solution = solution  # Create a pointer to the data

        # To do:  error handling

        return solution

    def _get_variable_name(self, ivar):
        variable_name_buffer = create_string_buffer(_VAR_BUFFER_SIZE)
        self._lib.GetVariableName(c_int(ivar), variable_name_buffer)

        return variable_name_buffer.value.decode("utf-8").strip()
