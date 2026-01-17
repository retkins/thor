""" Test sending and receiving immutable and mutable arrays to Rust
"""

from thor import count
import numpy as np 

def test_arrays():
    # Count the number of elements in a, return that as an int, and 
    # double each element, which is stored in b
    a = np.array([1,2,3.0])
    b = np.array([0.0, 0, 0])
    x = count(a,b)
    assert(x == 3) 
    assert((b - 2*a).all() == 0.0)