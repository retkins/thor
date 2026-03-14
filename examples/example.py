import thor
import numpy as np
from time import perf_counter
from matplotlib import pyplot as plt

size = 1.0
theta = 0.5
nthreads = 0
centroids, vol, jdensity = thor.testing.make_helmholtz(size)

# centroids = np.vstack((centroids, centroids))
# vol = np.hstack((vol, vol))
# jdensity = np.vstack((jdensity, jdensity))

print("Thor Example - Helmholtz Problem")
n = centroids.shape[0]
print(f"n = {n:.3e} ({n * n:.3e} interactions)")


start = perf_counter()
b = thor.bfield_octree(
    centroids, vol, jdensity, centroids, theta=theta, nthreads=nthreads
)
end = perf_counter()
elapsed = end - start

print(f"Elapsed time: {elapsed:.3f} s")
