""" Python bindings for thor
"""


from .biotsavart import bfield_direct, bfield_octree
from . import mesh
from . import test_utils
# import .test_utils as test_utils

__all__ = [
    "bfield_direct", 
    "bfield_octree", 
    "test_utils", 
    "mesh"
]