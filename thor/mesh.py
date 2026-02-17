""" Mesh generation and processing routines
"""

import numpy as np
from numpy.typing import NDArray 
from numpy import float64

def plot_mesh(x,y,z): 
    """ Make a scatter plot of element centroids
    """

    try: 
        import matplotlib.pyplot as plt 
        fig =plt.figure() 
        ax = fig.add_subplot(projection='3d')
        ax.scatter(x,y,z)
        plt.show()
    except ImportError:
        print("Error - matplotlib is not installed. Could not plot mesh.")


def mesh_step(infile: str, outfile: str, min_size: float, max_size: float): 
    """ Mesh a step file using gmsh
    """

    mshfile = infile.split('.')[0]+".msh"

    try:
        import gmsh
        gmsh.initialize()
        gmsh.model.occ.importShapes(infile)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", min_size)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", max_size)
        gmsh.model.mesh.generate(3)     # mesh 3d elements
        gmsh.write(mshfile)
        gmsh.finalize()
        print(f"Wrote gmsh mesh to `{mshfile}")
        process_elements(mshfile, outfile)

    except ImportError:
        print(f"Error - gmsh is not installed. Could not mesh file `{infile}`")


def mesh_step_tets(step_file: str, min_size: float, max_size: float, scale: float=1e-3) -> tuple[NDArray[float64], NDArray[float64], NDArray[float64]]:
    """
    Mesh a step file with gmsh and return tet element data. This is meant to 
    be used with the tet element source functionality.
    
    Args
    ---
    - `step_file`: Path to STEP file
    - `min_size`, `max_size`: Mesh element size bounds (in model units, usually mm)
    
    Returns
    ---
    (`nodes`, `centroids`, `volume`)
    - `nodes`: N*12-length flat array of nodal coordinates for each tet, row-major:
        [x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3, ...]
    - `centroids`: Nx3 array of the centroids of each element
    - `volume`: N-length array of volume of each element
    """

    import gmsh 

    # Setup gmsh and generate elements
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)  # suppress output
    gmsh.model.add("model")
    gmsh.model.occ.importShapes(step_file)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMin", min_size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", max_size)
    gmsh.model.mesh.generate(3)
    
    # Get all node coordinates: node_tags is 1-indexed
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    # node_coords is flat [x1,y1,z1, x2,y2,z2, ...]
    # Build a lookup from tag -> coordinates
    all_coords = node_coords.reshape(-1, 3)
    # node_tags might not be contiguous, so use a dict
    tag_to_idx = {int(tag): i for i, tag in enumerate(node_tags)}
    
    # Get tet elements (type 4 = linear tet with 4 nodes)
    tet_type = 4
    tet_tags, tet_node_tags = gmsh.model.mesh.getElementsByType(tet_type)
    
    n_tets = len(tet_tags)
    tet_connectivity = tet_node_tags.reshape(n_tets, 4)  # each row: 4 node tags
    
    # Build the 12*N flat node coordinate array
    nodes = np.zeros((n_tets, 4, 3))
    for i in range(n_tets):
        for j in range(4):
            idx = tag_to_idx[int(tet_connectivity[i, j])]
            nodes[i, j, :] = all_coords[idx]

    nodes *= scale
    
    # Centroids: average of 4 nodes (TODO: should this be done differently?)
    centroids = nodes.mean(axis=1) 
    
    # Volumes: V = |det([v1-v0, v2-v0, v3-v0])| / 6 like below, but now vectorized 
    v0 = nodes[:, 0, :]
    v1 = nodes[:, 1, :]
    v2 = nodes[:, 2, :]
    v3 = nodes[:, 3, :]
    
    d1 = v1 - v0
    d2 = v2 - v0
    d3 = v3 - v0
    
    cross = np.cross(d1, d2)
    det = np.sum(cross * d3, axis=1)
    volumes = np.abs(det) / 6.0
    
    # Flatten nodes to row-major 12*N
    nodes_flat = nodes.reshape(-1)  # [x0,y0,z0,x1,y1,z1,...] per tet
    
    gmsh.finalize()
    
    return nodes_flat, centroids, volumes


def tet_volume(p0, p1, p2, p3):
    """Calculate volume of a tetrahedron given 4 vertex coordinates:
        volume = |det(p1-p0, p2-p0, p3-p0)| / 6
    
        Each coordinate is a 3-length numpy array
        TODO: turn this into numpy operations for efficiency
    """

    v1 = p1 - p0
    v2 = p2 - p0
    v3 = p3 - p0
    return abs(np.dot(v1, np.cross(v2, v3))) / 6.0


def process_elements(infile: str, outfile: str, scale: float=1e-3):
    """ Convert gmsh .msh file into format readable by `thor` 
        and calculate the volume of each element
    """

    try:
        import gmsh
        gmsh.initialize()
        gmsh.open(infile)
        element_types, element_tags, node_tags = gmsh.model.mesh.getElements(3)

        with open(outfile, "w") as f:
            f.write("x,y,z,volume\n")
            
            for elem_type, elem_tags, elem_nodes in zip(element_types, element_tags, node_tags):

                elem_name, dim, order, num_nodes, local_coords, num_primary_nodes = gmsh.model.mesh.getElementProperties(elem_type)
                
                # Reshape node tags to have one row per element
                elem_nodes = elem_nodes.reshape(-1, num_nodes)
                
                print(f"Processing {len(elem_tags)} {elem_name} elements...")
                
                # Process each element
                for elem_tag, nodes in zip(elem_tags, elem_nodes):
                    # Get coordinates of all nodes in this element
                    coords = []
                    for node in nodes:
                        coord = gmsh.model.mesh.getNode(node)[0]
                        coords.append(coord)
                    coords = np.array(coords) * scale
                    
                    # Calculate centroid (average of node coordinates)
                    # TODO: is this the right way to calculate centroid?
                    centroid = np.mean(coords, axis=0)
                    
                    # Calculate volume based on element type
                    if "Tetrahedron" in elem_name:
                        volume = tet_volume(coords[0], coords[1], coords[2], coords[3])
                    
                    else:
                        print(f"Warning: Unknown element type {elem_name}, approximating volume")
                        volume = 0.0
                    
                    f.write(f"{centroid[0]:.6f},{centroid[1]:.6f},{centroid[2]:.6f},{volume:.10e}\n")

        print(f"Element data written to: {outfile}")
        gmsh.finalize()

    except ImportError:
        print(f"Error - gmsh is not installed. Could not process elements in file `{infile}`")
