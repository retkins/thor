""" Octree code in Julia for faster iteration
    
    Doesn't work yet...
"""

using LinearAlgebra, Printf

""" Biot Savart source point in 3D space
"""
struct Point 
    x::Float64; y::Float64; z::Float64 
    vJx::Float64; vJy::Float64; vJz::Float64 
    
    function Point(x::Real, y::Real, z::Real, vol::Real, Jx::Real, Jy::Real, Jz::Real)
        new(float(x), float(y), float(z), float(vol*Jx), float(vol*Jy), float(vol*Jz))
    end 
end


mutable struct Node 
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

"""
    find_octant(point::Point, node::Node)

Find the location in the Node where the Point should go 
Returns: integer value in the range [1,8] 
"""
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

"""
    calculate_span(points::Vector{Point})

Determine how big the volume is and where the centroid is 
Returns: (span, centroid), where
    span is the enclosing cube's side length 
    centroid is a 3-length vector of (xc, yc, zc)

TODO: should each octant be a cube? 
"""
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


# Add a leaf to the tree (node with a point source and no children)
function new_leaf!()
end

# 

"""
    build_octree(points::Vector{Point})

Build the octree from a collection of point sources
"""
function build_octree(points::Vector{Point})

    # Initial allocations for the nodes; assume nnodes=npts for now and reallocate later
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

    # Now we start at the root of the tree and descend to find a place for each source point
    for i in 1:npts 
        parent = 0 
        j_node = 1 
        point = points[i]

        while false
            node = nodes[j_node]

            # in which octant should the point go?
            octant = find_octant(point, node)
            j_child = node.children[octant]

            if j_child == 0 
                # There's no Node there currently, make a leaf Node

                # Establish a link between the current parent node and 
                # the leaf (child) node where the source point will be placed 
                node.children[octant] = current_child_pointer
                current_child_pointer += 1
                nodes[j_child].parent = j_node 

                # Now we've established the link, enter the new leaf node 
                # and update its values
                parent_node = node
                node = nodes[j_child]
                node.span = parent_node.span/2.0 


                
                nodes.ch

                # Update child pointer
                

                # We've placed the point, so now we can exit the loop
                break
            
            else
                # There's a node there already, but we need to do something different
                # if its a branch or a leaf 

                if nodes[j_child].point != 0 
                    # It's a leaf node, make it a branch and descend
                    # Now we need to:
                    #   1. turn a leaf into a branch
                    #   2. place the old leaf source point into a new leaf 
                    #   3. place the new leaf source point into a new leaf
                
                else 
                    # It's another branch node, descend further into the tree 
                    j_node = j_child 
                end

                
            end



        end
    end

    return Octree(points, nodes)

end