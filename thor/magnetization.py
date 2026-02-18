""" Operations for magnetic materials
"""
from typing import Callable

import numpy as np 
from numpy.typing import NDArray
from numpy import float64

from thor import LinearMaterial, MU0

def mag_force(
    centroids: NDArray[float64],
    vol: NDArray[float64],
    material: LinearMaterial,
    h_ext: Callable
) -> NDArray[float64]:
    """ Compute the force acting on the elements of a mesh with a magnetic material 
        under the influence of an external magnetic field. 

    Limitations:
    1. This does not (yet) consider the demagnetizing field 
    2. Materials are limited to linear (constant mu_r)

    Args:
        centroids: (N,3) centroid of each target element in 3d space
        vol: (N,) volume of each target element
        material: linear magnetic material
        h_ext(centroids): function to evaluate the external field 

    Returns:
        (N,3) force vector on each element
    """

    n: int = vol.shape[0]

    # keep the finite difference calculation within the element
    # perhaps make this constant across the mesh
    # delta = 0.5 * (3.0/(4.0*np.pi))*vol**(1.0/3.0)
    delta = 1e-4

    # grad_h[i,a,b] = dH_a/d_b
    grad_h = np.zeros((n, 3, 3)) 

    # second order finite differencing
    for axis in range(3):
        targets_plus = centroids.copy()
        targets_minus = centroids.copy()
        targets_plus[:, axis] += delta
        targets_minus[:, axis] -= delta
        
        # in the future, this should be:
        # h = h_ext + h_demag
        h_plus = h_ext(targets_plus)
        h_minus = h_ext(targets_minus)
        grad_h[:, :, axis] = (h_plus - h_minus) / (2.0 * delta)

    # in absence of demag field, assume moments are proportional to mu_r
    moments = material.chi(1) * vol[:,np.newaxis] * h_ext(centroids)

    # Force on each element
    forces = np.zeros((n, 3))
    for a in range(3):
        for b in range(3):
            forces[:, a] += MU0 * moments[:, b] * grad_h[:, b, a]   

    return forces