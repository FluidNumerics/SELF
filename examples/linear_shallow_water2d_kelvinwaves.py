from pyself.interface import SelfModel
from pyself.config import SelfModelConfig
import os
from datetime import datetime


pwd = os.path.dirname(os.path.abspath(__file__))
# Set the case directory to a new unique directory based in the time stamp
# get current working directory

case_directory = f"{os.getcwd()}/kelvinwaves-{datetime.now().strftime('%Y%m%d-%H%M%S')}"


def configure_geometry(config):
    # Configure geometry
    config.set_parameter(
        "geometry", "mesh_file", f"{pwd}/../share/mesh/Circle/Circle_mesh.h5"
    )
    config.set_parameter("geometry", "uniform_boundary_condition", "no_normal_flow")
    config.set_parameter("geometry", "control_degree", 7)
    config.set_parameter("geometry", "control_quadrature", "gauss")


def configure_time_options(config):
    # Configure time options
    config.set_parameter("time_options", "integrator", "euler")
    config.set_parameter("time_options", "dt", 0.0025)
    config.set_parameter("time_options", "start_time", 0.0)
    config.set_parameter("time_options", "duration", 1.0)
    config.set_parameter("time_options", "io_interval", 0.05)
    config.set_parameter("time_options", "update_interval", 50)


def configure_shallow_water(config):
    # Configure shallow water parameters
    config.set_parameter("linear-shallow-water-2d", "g", 1.0)
    config.set_parameter("linear-shallow-water-2d", "H", 1.0)
    config.set_parameter("linear-shallow-water-2d", "Cd", 0.25)
    config.set_parameter("linear-shallow-water-2d", "f0", 10.0)
    config.set_parameter("linear-shallow-water-2d", "beta", 0.0)


def main():

    config = SelfModelConfig(case_directory=case_directory)
    config.config["model_name"] = "linear-shallow-water-2d"

    configure_geometry(config)
    configure_time_options(config)
    configure_shallow_water(config)

    # Create the model
    model = SelfModel(config=config)

    x, y = model.get_coordinates()


if __name__ == "__main__":
    main()
