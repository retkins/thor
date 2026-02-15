""" Python bindings for thor
"""


from .biotsavart import bfield_direct, bfield_octree, bfield_dualtree, bfield_hexahedron
from . import mesh
from . import test_utils
# import .test_utils as test_utils

__all__ = [
    "bfield_direct", 
    "bfield_octree", 
    "test_utils", 
    "mesh", 
    "bfield_dualtree", 
    "bfield_hexahedron"
]