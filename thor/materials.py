""" Define material properties for magnetostatics calculations
"""


from dataclasses import dataclass
from numpy import pi
from numpy.typing import NDArray 
from numpy import float64, interp
import thor

""" Magnetic permeability of free space
"""
MU0: float = (4.0*pi) * 10**-7 

@dataclass
class BHCurve:
    """ B-H curve for a nonlinear magnetic material
    """ 
    h_values: NDArray[float64]
    b_values: NDArray[float64]

    def lookup(self, h: float) -> float:
        """ Linearly interpolate from the B-H curve
        """

        return interp(h, self.h_values, self.b_values)


class Material:

    def mu_r(self, h: float) -> float: 
        raise NotImplementedError 

    def chi(self, h: float) -> float: 
        """ Return the volume susceptibility of the material, defined by
            `chi = mu_r - 1.0`
        """
        return self.mu_r(h) - 1.0 


@dataclass
class FreeSpace(Material):

    def mu_r(self, h: float) -> float:
        return 1.0 


@dataclass
class LinearMaterial(Material):
    """ A linear magnetic material that has a constant mu_r for all
        values of applied H-field.
    """

    _mu_r: float 

    def __init__(self, _mu_r):
        self._mu_r = _mu_r 

    def mu_r(self, h: float): 
        return self._mu_r 


@dataclass
class NonlinearMaterial(Material):
    """ A nonlinear magnetic material that has a B-H curve
    """

    curve: BHCurve 

    def __init__(self, curve: BHCurve):
        self.curve = curve 

    def mu_r(self, h: float) -> float:
        b: float = self.curve.lookup(h)

        if h != 0.0:
            return (b/thor.MU0)/h - 1.0
        else: 
            return 1.0