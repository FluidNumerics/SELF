"""
self

A python package for interfacing with the Spectral Element Library in Fortran (https://github.com/fluidnumerics/SELF).
"""

from ._version import version

__version__ = version
__title__ = "pyself"
__author__ = "Dr. Joe Schoonover"
__credits__ = "Fluid Numerics LLC"


from pyself.config import *
from pyself.geometry import *
from pyself.interface import *
from pyself.lagrange import *
from pyself.model import *
