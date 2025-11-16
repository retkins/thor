""" Octree code in Julia for faster iteration
"""

using LinearAlgebra, Printf

struct Point 
    x::Float64
    y::Float64 
    z::Float64 
    vJx::Float64
    vJy::Float64 
    vJz::Float64 
    
    function Point(x::Real, y::Real, z::Real, vol::Real, Jx::Real, Jy::Real, Jz::Real)
        new(float(x), float(y), float(z), float(vol*Jx), float(vol*Jy), float(vol*Jz))
    end 
end


struct Node 
    span::Float64                                   # full width of a side
    x::Float64; y::Float64; z::Float64
    point::Int64
    parent::Int64
    children::Vector{Int64}
    vJx::Float64; vJy::Float64; vJz::Float64        # volume*J
    mx::Float64; my::Float64; mz::Float64           # magnetic moment
    function Node(span::Real, x::Real, y::Real, z::Real, parent::Integer)
        new(float(span), float(x), float(y), float(z), 0, parent, zeros(8), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
end

function initialize_nodes(nnodes)
    nodes = Vector{Node}(undef, nnodes)
    for i in 1:nnodes 
        nodes[i] = Node(0, 0, 0, 0, 0) 
    end
    return nodes
end

struct Octree
    points::Vector{Point}
    nodes::Vector{Node}
end

# figure out where the point should go in the node 
function find_octant(point::Point, node::Node)
    # convenience definitions 
    x = point.x; y = point.y; z = point.z 
    cx = node.x; cy = node.y; cz = node.z 

    if (x >= cx) && (y >= cy) && (z >= cz) 
        return 1 
    elseif (x < cx) && (y >= cy) && (z >= cz)
        return 2
    elseif (x < cx) && (y < cy) && (z >= cz)
        return 3 
    elseif (x >= cx) && (y < cy) && (z >= cz)
        return 4 
    elseif (x >= cx) && (y >= cy) && (z < cz) 
        return 1 
    elseif (x < cx) && (y >= cy) && (z < cz)
        return 2
    elseif (x < cx) && (y < cy) && (z < cz)
        return 3 
    elseif (x >= cx) && (y < cy) && (z < cz)
        return 4 
    end

end

# figure out how big the space is; TODO: should each octant be a cube?
function calculate_span(points::Vector{Point})
    xmax = 0.0; xmin = 0.0; ymax = 0.0; ymin = 0.0; zmax = 0.0; zmin = 0.0 
    x = point.x; y = point.y; z = point.z 

    for i in 1:length(points) 
        if x > xmax 
            xmax = x 
        elseif x < xmin 
            xmin = x 
        
        if y > ymax 
            ymax = y 
        elseif y < ymin 
            ymin = y 

        if z > zmax 
            zmax = z 
        elseif z < zmin 
            zmin = z 
    end

    max_size = maximum([xmax, ymax, zmax])
    min_size = minimum([xmin, ymin, zmin])
    span = max_size - min_size
    center = 0.5 * (max_size + min_size)    

    return span, [center, center, center]
end

# Build the octree
function build_octree(points::Vector{Point})

    npts = length(points) 
    nnodes = npts
    nodes = initialize_nodes(nnodes)       # may need to be resized
    
    # Get info about the points and make root node (first child = 2)
    span, centroid = calculate_span(points)
    parent = 0 
    j_node = 1
    current_child_pointer = 2
    nodes[j_node] = Node(span, centroid[1], centroid[2], centroid[3], parent)
    node = nodes[j_nodes]
    for i in 1:npts 

        while false

            # in which octant should the point go?
            octant = find_octant(points[i], node)

            # Check if there's a Node there 
            if nodes.children[octant] == 0 
                # No node, make one and add the point to it (leaf)
                nodes.children[octant] = current_child_pointer

                # Update child pointer
                current_child_pointer += 1

                break
            
            # There's a node there with a point (no children)
            if nodes[j_octant]. == 0 



        end
    end

    return Octree(points, nodes)

end