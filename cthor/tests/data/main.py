import numpy as np
import matplotlib.pyplot as plt 
import gmsh

filename: str = "solenoid.stp" 
max_size: float = 3.5 
min_size: float = 1.0
convert_to_meters: bool = True

def plot(x,y,z): 
    fig =plt.figure() 
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x,y,z)
    plt.show()


def mesh_step(file: str, min_size: float, max_size: float): 
    """ Mesh a step file using gmsh
    """
    gmsh.initialize()
    gmsh.model.occ.importShapes(filename)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", min_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", max_size)
    gmsh.model.mesh.generate(3)     # mesh 3d elements
    gmsh.write(file+".msh")
    gmsh.finalize()

def tet_volume(p0, p1, p2, p3):
    """Calculate volume of a tetrahedron given 4 vertex coordinates:
        volume = |det(p1-p0, p2-p0, p3-p0)| / 6
    
        Each coordinate is a 3-length numpy array
    """

    v1 = p1 - p0
    v2 = p2 - p0
    v3 = p3 - p0
    return abs(np.dot(v1, np.cross(v2, v3))) / 6.0


def process_elements(file: str):

    gmsh.initialize()
    gmsh.open(file+".msh")
    element_types, element_tags, node_tags = gmsh.model.mesh.getElements(3)

    with open(file+"_data.csv", "w") as f:
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
                coords = np.array(coords)
                
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

    print("Element data written to: "+file+"_data.csv")
    gmsh.finalize()


def write_csv(file: str):
    data = np.loadtxt(file+"_data.csv", delimiter=',', skiprows=1)
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    vol = data[:,3] 
    if convert_to_meters: 
        x *= 1e-3
        y *= 1e-3
        z *= 1e-3
        vol *= 1e-9
    print(f"Max z = {np.max(z)}")
    # print(f"Total volume: {np.sum(vol)} mm3")

    mu0 = 4*np.pi*10**-7 
    b = 5 
    # L = 0.25 
    thk = 25e-3
    J = b /(mu0 * thk)
    print(f"J = {J:.3e}")

    theta = np.atan2(y,x) 
    jx = -J*np.sin(theta) 
    jy = J*np.cos(theta) 
    jz = 0.0*theta

    output_data = np.vstack((x,y,z,vol,jx,jy,jz)).T 
    np.savetxt(file+".csv", output_data, delimiter=',')


def main():

    file = filename.split(".")[0]
    mesh_step(file, min_size, max_size)
    process_elements(file) 
    write_csv(file)



if __name__ == "__main__":
    main()
