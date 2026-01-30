""" Mesh generation and processing routines
"""

import numpy as np

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
