""" Python bindings for thor
"""


from ._thor import test, count_rs, _bfield_direct, _bfield_octree
import numpy as np
from numpy import float64, ascontiguousarray, zeros, hstack, newaxis
from numpy.typing import NDArray

Nx3Arrray = NDArray[float64]


def test2(x): 
    """ Confirm that we can call Python functions without name collisions
    """

    return x - 1.0


def gemv(a: NDArray[float64], b: NDArray[float64], alpha: float, beta: float, y: NDArray[float64]) -> NDArray[float64]:
    """ Calculate matrix-vector product in Rust
    
    This function confirms both 1d and 2d, immutable and mutable numpy
    matrices work properly
    """ 

    _thor.gemv(a, b, alpha, beta, y)
    return y

def count(a, b) -> int: 
    return count_rs(a, b)


def bfield_direct(
    centroids: NDArray[float64], 
    vol: NDArray[float64], 
    jdensity: NDArray[float64], 
    targets: NDArray[float64], 
) -> NDArray[float64]: 

    ntargets = targets.shape[0]
    bx = zeros(ntargets)
    by = zeros(ntargets)
    bz = zeros(ntargets)


    _bfield_direct(
        ascontiguousarray(centroids[:,0]),
        ascontiguousarray(centroids[:,1]),
        ascontiguousarray(centroids[:,2]),
        ascontiguousarray(vol[:]),
        ascontiguousarray(jdensity[:,0]),
        ascontiguousarray(jdensity[:,1]),
        ascontiguousarray(jdensity[:,2]),
        ascontiguousarray(targets[:,0]),
        ascontiguousarray(targets[:,1]),
        ascontiguousarray(targets[:,2]),
        ascontiguousarray(bx[:]),
        ascontiguousarray(by[:]),
        ascontiguousarray(bz[:])
    )

    return hstack((bx[:, newaxis], by[:, newaxis], bz[:, newaxis]))


def bfield_octree(
    centroids: NDArray[float64], 
    vol: NDArray[float64], 
    jdensity: NDArray[float64], 
    targets: NDArray[float64], 
    theta: float=0.5, 
    leaf_threshold: int=1
) -> NDArray[float64]: 

    ntargets = targets.shape[0]
    bx = zeros(ntargets)
    by = zeros(ntargets)
    bz = zeros(ntargets)

    _bfield_octree(
        ascontiguousarray(centroids[:,0]),
        ascontiguousarray(centroids[:,1]),
        ascontiguousarray(centroids[:,2]),
        ascontiguousarray(vol[:]),
        ascontiguousarray(jdensity[:,0]),
        ascontiguousarray(jdensity[:,1]),
        ascontiguousarray(jdensity[:,2]),
        ascontiguousarray(targets[:,0]),
        ascontiguousarray(targets[:,1]),
        ascontiguousarray(targets[:,2]),
        ascontiguousarray(bx[:]),
        ascontiguousarray(by[:]),
        ascontiguousarray(bz[:]), 
        theta, 
        leaf_threshold
    )

    return hstack((bx[:, newaxis], by[:, newaxis], bz[:, newaxis]))

__all__ = [
    "test", 
    "test2", 
    "gemv"
]