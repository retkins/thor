""" Python bindings for thor
"""

from ._thor import test

def test2(x): 
    """ Confirm that we can call Python functions without name collisions
    """

    return x - 1.0


__all__ = [
    "test", 
    "test2"
]