import json
from typing import Optional, Dict, Any

class SelfModelConfig:
    def __init__(self, config_file: Optional[str] = None, case_directory: Optional[str] = None):
        """Initialize the SELF model configuration from a JSON file or defaults."""
        self.config = self.default_config()

        self.config_file = config_file
        
        if config_file:
            self.load_config(config_file)

        if case_directory:
            self.case_directory = case_directory
        else:
            self.case_directory = "."

        

    @staticmethod
    def default_config() -> Dict[str, Any]:
        """Return default configuration based on the JSON schema."""
        return {
            "version": "v0.0.0",
            "model_name": "linear-shallow-water-2d",
            "geometry": {
                "mesh_file": "",
                "uniform_boundary_condition": "no_normal_flow",
                "control_degree": 7,
                "control_quadrature": "gauss",
                "target_degree": 10,
                "target_quadrature": "uniform",
                "nX": 5,
                "nY": 5,
                "nZ": 5,
                "nTx": 1,
                "nTy": 1,
                "nTz": 1,
                "dx": 0.02,
                "dy": 0.02,
                "dz": 0.02
            },
            "time_options": {
                "integrator": "euler",
                "dt": 0.001,
                "cfl_max": 0.5,
                "start_time": 0.0,
                "duration": 1.0,
                "io_interval": 0.1,
                "update_interval": 50
            },
            "units": {
                "time": "s",
                "length": "m",
                "mass": "kg"
            },
            "linear-shallow-water-2d": {
                "g": 1.0,
                "H": 1.0,
                "Cd": 0.01,
                "f0": 0.0,
                "beta": 0.0,
                "initial_conditions": {
                    "geostrophic_balance": false,
                    "file": ""
                    "u": 0.0,
                    "v": 0.0,
                    "eta": 0.0
                },
                "boundary_conditions": {
                    "time_deppendent": false,
                    "dt": 0.0,
                    "from_initial_conditions": false,
                    "u": 0.0,
                    "v": 0.0,
                    "eta": 0.0
                }
            }
        }

    def load_config(self, file_path: str):
        """Load configuration from a JSON file."""
        with open(file_path, "r") as f:
            self.config.update(json.load(f))

    def save_config(self, file_path: str):
        """Save configuration to a JSON file."""
        with open(file_path, "w") as f:
            json.dump(self.config, f, indent=4)

    def set_parameter(self, section: str, key: str, value: Any):
        """Set a specific parameter within the configuration."""
        if section in self.config and key in self.config[section]:
            self.config[section][key] = value
        else:
            raise KeyError(f"Invalid section '{section}' or key '{key}'.")

    def get_parameter(self, section: str, key: str) -> Any:
        """Retrieve a specific parameter value."""
        return self.config.get(section, {}).get(key, None)


# Example Usage
config = SelfModelConfig()
config.set_parameter("geometry", "nX", 10)
config.set_parameter("time_options", "dt", 0.005)
print(config.generate_fortran_input())
