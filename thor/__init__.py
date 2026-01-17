""" Python bindings for thor
"""


from ._thor import test
from ._thor import count_rs
import numpy as np
from numpy import float64
from numpy.typing import NDArray


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

__all__ = [
    "test", 
    "test2", 
    "gemv"
]