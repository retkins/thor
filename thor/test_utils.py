""" Utility functions for tests
"""

import numpy as np 
from numpy import array, float64 
from numpy.typing import NDArray

def mean_squared_error(baseline: NDArray[float64], measurement: NDArray[float64]) -> float:
    """ Compute the mean squared error of `measurement` measured against `baseline`

    Error is computed according to this reference:
    https://en.wikipedia.org/wiki/Mean_squared_error#Predictor

    This is an absolute, not a relative error measurement.

    Args:
    ---
        baseline: N-length array of 'ground-truth' values 
        measurement: N-length array of values to which determine error (deviation) from baseline

    """

    assert(baseline.shape == measurement.shape)
    assert(len(baseline.shape) == 1)
    n: int = baseline.shape[0]
    assert(n > 0)

    return (1/n) * np.sum((baseline - measurement)**2)


def mean_relative_error(baseline: NDArray[float64], measurement: NDArray[float64]) -> float:
    """ Compute the mean *relative* error of `measurement` against `baseline`
    """

    assert(baseline.shape == measurement.shape)
    assert(len(baseline.shape) == 1)
    n: int = baseline.shape[0]
    assert(n > 0)

    relative_diff = (measurement - baseline) / baseline 
    return (1/n) * np.sum(np.abs(relative_diff)) 


def mean_absolute_error(baseline: NDArray[float64], measurement: NDArray[float64]) -> float:
    """ Compute the mean absolute error of `measurement` against `baseline`
    """

    assert(baseline.shape == measurement.shape)
    assert(len(baseline.shape) == 1)
    n: int = baseline.shape[0]
    assert(n > 0)

    return (1/n) * np.sum(baseline - measurement)


def smape(baseline: NDArray[float64], measurement: NDArray[float64]) -> float:
    """ Compute the symmetric mean absolute percentage error of `measurement` against `baseline`

    SMAPE is defined here:
    https://en.wikipedia.org/wiki/Symmetric_mean_absolute_percentage_error
    """

    assert(baseline.shape == measurement.shape)
    assert(len(baseline.shape) == 1)
    n: int = baseline.shape[0]
    assert(n > 0)

    numerator = np.abs(measurement - baseline) 
    denominator = np.abs(measurement) + np.abs(baseline)
    return (2/n) * np.sum(numerator / denominator)

