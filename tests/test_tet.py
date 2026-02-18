""" 
"""

import thor 
import numpy as np 
import matplotlib.pyplot as plt
from time import perf_counter

min_size: float = 2.5 
max_size: float = 10.0
jmag: float = 1e8 
theta_pt: float = 0.25
theta_tet: float = 0.25
nthreads: int = 1
ntargets_axis: int = 100

# load mesh
nodes, centroids, vol = thor.mesh.mesh_step_tets("tests/data/solenoid.stp", min_size, max_size)
n = len(vol)

print("Thor: Solenoid Test for tetrahedron sources in octree \n---")
print(f"\tn = {n}")
print(f"\ttheta (pt) =  {theta_pt:.3f}")
print(f"\ttheta (tet) = {theta_tet:.3f}")
print(f"\tNthreads = {nthreads}")
print("\tErrors reported relative to direct pt source summation\n")

# assign current density
jdensity = np.zeros((n,3))
phi = np.atan2(centroids[:,1], centroids[:,0])
jdensity[:,0] = -jmag*np.sin(phi)
jdensity[:,1] = jmag*np.cos(phi)


# --- 
# Solution at center of solenoid
# ---

targets_axis = np.zeros((ntargets_axis, 3))
targets_axis[:,2] = np.linspace(-0.125, 0.125, ntargets_axis)
# targets_axis[:,0] = np.linspace(0, 0.10, ntargets_axis)

bdirect_pt_axis = thor.bfield_direct(centroids, vol, jdensity, targets_axis)
bdirect_tet_axis = thor.bfield_tetrahedrons_direct(nodes, centroids, vol, jdensity, targets_axis, nthreads=nthreads)
boctree_tet_axis = thor.bfield_tetrahedrons(nodes, centroids, vol, jdensity, targets_axis, theta=theta_tet, nthreads=nthreads)
boctree_pt_axis = thor.bfield_octree(centroids, vol, jdensity, targets_axis, nthreads=nthreads,theta=theta_pt)

# Errors along axis
print("Bfield at solenoid center: ")
i: int = ntargets_axis//2
print(f"\tDirect (pt) solution:  ({bdirect_pt_axis[i,0]:.6f}, {bdirect_pt_axis[i,1]:.6f}, {bdirect_pt_axis[i,2]:.6f})")
print(f"\tDirect (tet) solution: ({bdirect_tet_axis[i,0]:.6f}, {bdirect_tet_axis[i,1]:.6f}, {bdirect_tet_axis[i,2]:.6f})")
print(f"\tOctree (pt) solution:  ({boctree_pt_axis[i,0]:.6f}, {boctree_pt_axis[i,1]:.6f}, {boctree_pt_axis[i,2]:.6f})")
print(f"\tOctree (tet) solution: ({boctree_tet_axis[i,0]:.6f}, {boctree_tet_axis[i,1]:.6f}, {boctree_tet_axis[i,2]:.6f})")
err = (bdirect_pt_axis[i,2] - boctree_pt_axis[i,2])/bdirect_pt_axis[i,2]
print(f"Error at center: {100*err:.3f} %")

bmag_direct_axis = np.linalg.norm(bdirect_pt_axis, axis=1) 
bmag_octree_axis = np.linalg.norm(boctree_pt_axis, axis=1)
err_axis = thor.test_utils.smape(bmag_direct_axis, bmag_octree_axis)
print(f"Mean error along solenoid axis (|z| < radius): {err_axis*100:.2}%")

#
# Solve for self-fields
# 

targets = centroids
ntargets = targets.shape[0]
nsources = n

start = perf_counter()
bdirect_pt = boctree = thor.bfield_direct(centroids, vol, jdensity, targets, nthreads=nthreads)
end = perf_counter() 
direct_pt_elapsed = end - start

start = perf_counter()
bdirect_tet = thor.bfield_tetrahedrons_direct(nodes, centroids, vol, jdensity, targets)
end = perf_counter() 
direct_tet_elapsed = end - start

start = perf_counter()
boctree_pt = thor.bfield_octree(centroids, vol, jdensity, targets, nthreads=nthreads,theta=theta_pt)
end = perf_counter() 
octree_pt_elapsed = end - start 

start = perf_counter()
boctree_tet = thor.bfield_tetrahedrons(nodes, centroids, vol, jdensity, targets, theta=theta_tet, nthreads=nthreads)
end = perf_counter() 
octree_tet_elapsed = end - start 

# Errors on mesh
bmag_direct_pt = np.linalg.norm(bdirect_pt, axis=1) 
bmag_direct_et = np.linalg.norm(bdirect_tet, axis=1) 
bmag_octree_pt = np.linalg.norm(boctree_pt, axis=1) 
bmag_octree_tet = np.linalg.norm(boctree_tet, axis=1) 
err_mesh_pt = thor.test_utils.smape(bmag_direct_pt, bmag_octree_pt)
err_mesh_tet = thor.test_utils.smape(bmag_direct_pt, bmag_octree_tet)
print(f"Mean fields error within the mesh (pt sources):  {err_mesh_pt*100:.2}%")
print(f"Mean fields error within the mesh (tet sources): {err_mesh_tet*100:.2}%")

print("Times: ")
print(f"\tDirect (pt) solution time:  {direct_pt_elapsed*1e3:.3f} ms")
print(f"\tDirect (tet) solution time: {direct_tet_elapsed*1e3:.3f} ms")
print(f"\t\tSpeedup: {direct_pt_elapsed/direct_tet_elapsed:.2f}x")
print(f"\tOctree (pt) solution time:  {octree_pt_elapsed*1e3:.3f} ms")
print(f"\t\tSpeedup: {direct_pt_elapsed/octree_pt_elapsed:.2f}x")
print(f"\tOctree (tet) solution time: {octree_tet_elapsed*1e3:.3f} ms")
print(f"\t\tSpeedup: {direct_pt_elapsed/octree_tet_elapsed:.2f}x")


fig = plt.figure()
ax = fig.add_subplot()
ax.plot(targets_axis[:,2], bdirect_pt_axis[:,2], 'k', label="Direct (pt)")
ax.plot(targets_axis[:,2], bdirect_tet_axis[:,2], 'b--', label="Direct (tet)")

ax.plot(targets_axis[:,2], boctree_pt_axis[:,2], 'r.', label="Octree (pt)")
ax.plot(targets_axis[:,2], boctree_tet_axis[:,2], 'g^', label="Octree (tet)")

ax.set_xlabel("Distance from Solenoid Center (Z-axis) [m]")
ax.set_ylabel("Field Along Solenoid Axis (Bz) [T]")
ax.set_title("Thor - Solenoid Test - Pt vs Tet Sources")
plt.legend()
plt.savefig("tests/fig/solenoid_tet_test.png")
# print(boctree_tet_axis)