import numpy as np 
import matplotlib.pyplot as plt
import thor 
from time import perf_counter

data = np.loadtxt("tests/solenoid_72255.csv", delimiter=',')

centroids = data[:,0:3]
vol = data[:,3]
jdensity = data[:, 4:7]

nsources = centroids.shape[0]
ntargets = nsources

targets = np.zeros((ntargets, 3))
targets[:,2] = np.linspace(0.0, 0.5, ntargets)

# targets = centroids
# ntargets = targets.shape[0]

start = perf_counter()
bdirect = boctree = thor.bfield_direct(centroids, vol, jdensity, targets)
end = perf_counter() 
direct_elapsed = end - start

start = perf_counter()
boctree = thor.bfield_octree(centroids, vol, jdensity, targets, theta=0.25, leaf_threshold=1)
end = perf_counter() 
octree_elapsed = end - start 

print("Thor: mesh solenoid test")
print(f"Problem size: {nsources} x {ntargets} ({nsources*ntargets:.3e} interactions)")
print("Bfield at solenoid center: ")
print(f"Direct solution: ({bdirect[0,0]:.6f}, {bdirect[0,1]:.6f}, {bdirect[0,2]:.6f})")
print(f"Octree solution: ({boctree[0,0]:.6f}, {boctree[0,1]:.6f}, {boctree[0,2]:.6f})")
err = (bdirect[0,2] - boctree[0,2])/bdirect[0,2]
print(f"Error: {100*err:.3f} %")

print("Times: ")
print(f"Direct solution time: {direct_elapsed*1e3:.3f} ms")
print(f"Octree solution time: {octree_elapsed*1e3:.3f} ms")
print(f"Speedup: {direct_elapsed/octree_elapsed:.2f}x")

# for theta in [0.1, 0.3, 0.5, 0.7, 0.9]:
#     start = perf_counter()
#     boctree = thor.bfield_octree(centroids, vol, jdensity, targets, theta=theta, leaf_threshold=32)
#     elapsed = perf_counter() - start
#     err = (bdirect[0,2] - boctree[0,2])/bdirect[0,2]
#     print(f"theta={theta}: error={100*err:.3f}%, time={elapsed*1e3:.1f}ms")

fig = plt.figure()
ax = fig.add_subplot()
ax.plot(targets[:,2], bdirect[:,2], label="Direct")
ax.plot(targets[:,2], boctree[:,2], 'rs', label="Octree")
ax.set_xlabel("Distance from Solenoid Center (Z-axis) [m]")
ax.set_ylabel("Field Along Solenoid Axis (Bz) [T]")
ax.set_title("Thor - Mesh Solenoid Test")
# plt.show()