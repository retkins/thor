import gmsh
import numpy as np
import thor
import time

infile: str = "tests/data/sphere.stp"
outfile: str = "tests/data/sphere.data"
min_size: float = 10.0
max_size: float = 15.0
thor.mesh.mesh_step(infile, outfile, min_size, max_size)

gmsh.initialize()
gmsh.open("tests/data/sphere.msh")

# Get all nodes: returns (tags, coords, parametricCoords)
node_tags, coords, _ = gmsh.model.mesh.getNodes()

# coords is flat [x0,y0,z0,x1,y1,z1,...], reshape to (Nn, 3)
nodes = np.array(coords).reshape(-1, 3)

# Build compact renumbering: gmsh tags can be sparse/non-sequential
tag_to_compact = {tag: i for i, tag in enumerate(node_tags)}

# Get tet elements (type 4 = 4-node tetrahedra)
elem_tags, elem_node_tags = gmsh.model.mesh.getElementsByType(4)

# elem_node_tags is flat [n0,n1,n2,n3, n0,n1,n2,n3, ...], reshape to (Ne, 4)
raw_connectivity = np.array(elem_node_tags).reshape(-1, 4)

# Renumber to compact 0-based indices
connectivity = np.array([[tag_to_compact[tag] for tag in elem] for elem in raw_connectivity], dtype=np.uint32)
gmsh.finalize()

mat = thor.materials.LinearMaterial(1.4)
h_external = np.zeros((connectivity.shape[0], 3))
b_ext_mag: float = 1.0
h_ext_mag: float = b_ext_mag / thor.MU0
h_external[:, 2] = b_ext_mag / thor.MU0

start = time.perf_counter()
M, Htotal = thor.magnetization.demag_tet4(nodes, connectivity, mat, h_external, nthreads_requested=6)
elapsed = time.perf_counter() - start

Mnorm = np.linalg.norm(M, axis=1)
Hnorm = np.linalg.norm(Htotal, axis=1)

print("Demag tet test - magnetized sphere")
print(f"{nodes.shape[0]} nodes, {connectivity.shape[0]} elements")
print(f"Total solver time: {elapsed:.3f} sec")
print(f"Avg/min/max M: {np.average(Mnorm)} / {np.min(Mnorm):.3e} / {np.max(Mnorm):.3e}")
print(f"Avg/min/max H: {np.average(Hnorm)} / {np.min(Hnorm):.3e} / {np.max(Hnorm):.3e}")

h_z_analytical = h_ext_mag / (1 + mat.chi(1) / 3)
h_z_mean = np.average(Htotal[:, 2])


def test_magnetized_sphere():
    err = float(np.abs(h_z_mean - h_z_analytical) / h_z_analytical)
    print(f"Error = {err * 100:.3f} %")
    assert err < 0.01


if __name__ == "__main__":
    test_magnetized_sphere()
    print("Test passed!")
