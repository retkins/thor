""" Python bindings for thor
"""


from .biotsavart import (
    bfield_direct, 
    bfield_octree, 
    bfield_dualtree, 
    bfield_hexahedron, 
    bfield_tetrahedrons, 
    bfield_tetrahedrons_direct, 
    hfield_dipole, 
    hfield_dipole_tetrahedrons
)
from .materials import MU0, FreeSpace, LinearMaterial, NonlinearMaterial, BHCurve
from . import mesh
from . import test_utils
from .magnetization import mag_force


__all__ = [
    "MU0",
    "FreeSpace", 
    "LinearMaterial", 
    "NonlinearMaterial",
    "BHCurve",
    "bfield_direct", 
    "bfield_octree", 
    "test_utils", 
    "mesh", 
    "bfield_dualtree", 
    "bfield_hexahedron",
    "bfield_tetrahedrons",
    "bfield_tetrahedrons_direct",
    "hfield_dipole",
    "hfield_dipole_tetrahedrons",
    "mag_force"
]