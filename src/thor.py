""" Python bindings for thor 
"""

from enum import Enum
import numpy as np 
import matplotlib.pyplot as plt
from ctypes import CDLL, c_int, c_double, c_size_t
import sys, os
from numpy.ctypeslib import ndpointer
from numpy import ascontiguousarray
NDArray = np.ndarray[np.float64]
darray = ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")

this_file = os.path.abspath(__file__) 
build_dir = os.path.join(os.path.dirname(this_file), os.path.pardir, "build")
files = os.listdir(build_dir)

class Method(Enum):
    Direct = 1
    Octree = 2 

try: 
    if "thorlib.so" not in files:
        raise FileNotFoundError
    
    clib_file = os.path.join(build_dir, "thorlib.so")
    clib = CDLL(clib_file)
    _bfield_direct = clib.bfield_direct_simd
    _bfield_direct.argtypes = [darray, darray, darray, darray, darray, darray, darray, c_size_t, darray, darray, darray, c_size_t, darray, darray, darray, c_int]
    _bfield_direct.restype = c_int
    _bfield_octree = clib.bfield_octree
    _bfield_octree.argtypes = [darray, darray, darray, darray, darray, darray, darray, c_size_t, darray, darray, darray, c_size_t, darray, darray, darray, c_int, c_double]
    _bfield_octree.restype = c_int

except FileNotFoundError:
    print("Error! Could not find thor lib. Perhaps try building the program?")

# print(files)

def bfield(targets:NDArray, sources: NDArray, method: Method=Method.Direct, nthreads: int=1, phi: float=0.5) -> NDArray: 

    n_sources = sources.shape[0] 
    n_targets = targets.shape[0]
    Bx = ascontiguousarray(np.zeros(n_targets, dtype=np.float64))
    By = ascontiguousarray(np.zeros(n_targets, dtype=np.float64))
    Bz = ascontiguousarray(np.zeros(n_targets, dtype=np.float64))

    if method == Method.Direct:
        success = _bfield_direct(
            ascontiguousarray(sources[:,0]), 
            ascontiguousarray(sources[:,1]), 
            ascontiguousarray(sources[:,2]), 
            ascontiguousarray(sources[:,3]), 
            ascontiguousarray(sources[:,4]), 
            ascontiguousarray(sources[:,5]), 
            ascontiguousarray(sources[:,6]), 
            c_size_t(n_sources), 
            ascontiguousarray(targets[:,0]), 
            ascontiguousarray(targets[:,1]), 
            ascontiguousarray(targets[:,2]), 
            c_size_t(n_targets), 
            Bx, By, Bz,
            c_int(nthreads)
        )
    else:
        success = _bfield_octree(
            ascontiguousarray(sources[:,0]), 
            ascontiguousarray(sources[:,1]), 
            ascontiguousarray(sources[:,2]), 
            ascontiguousarray(sources[:,3]), 
            ascontiguousarray(sources[:,4]), 
            ascontiguousarray(sources[:,5]), 
            ascontiguousarray(sources[:,6]), 
            c_size_t(n_sources), 
            ascontiguousarray(targets[:,0]), 
            ascontiguousarray(targets[:,1]), 
            ascontiguousarray(targets[:,2]), 
            c_size_t(n_targets), 
            Bx, By, Bz,
            c_int(nthreads), c_double(phi)
        )


    return np.hstack((Bx[:,np.newaxis], By[:,np.newaxis], Bz[:,np.newaxis]))

data = np.loadtxt("tests/data/solenoid.csv", delimiter=',')
targets = np.zeros((10,3))
targets[:,0] = 0.0
sources = data
bdirect = bfield(targets, sources, method=Method.Direct)
boctree = bfield(targets, sources, method=Method.Octree, phi=1e12)
print(f"Error = {100*(boctree[0,2] - bdirect[0,2])/bdirect[0,2]:.2f}%")
print(f"Bz direct = {bdirect[0,2]:.4f}, Bz octree = {boctree[0,2]:.4f}")