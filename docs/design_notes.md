# Design Notes for Thor

## Tests and Benchmarks
Define the following tests for accuracy checks and performance benchmarks:
1. Random point cloud in 3D space: compare the self-field at each point to the self-field computed by direct summation. This will result in a balanced octree.
2. Ring: create a circular ring around the origin for the source points. Output the result on a line along the central axis. This compares against a theoretical solution that does not involve self-field.
3. Random point cloud in 3D space; unbalanced octree. Randomly assign points in an uneven fashion throughout the top-level octants. Determine speed benefits of the octree code when the tree is not balanced.

## Limitations
Indexes are performed using unsigned 32-bit integers; therefore, the maximum number of nodes and source/target points in a given tree is 2^32 (~4.3B). 
For trees that must contain that many nodes, please use another software code.

## Features to add in 'refactor' branch
### Immediate
- [ ] Replace pointers with array indices
- [ ] Replace recursion with explicit traversal
- [ ] Add parallel processing for direct and octree methods (this will also determine if stack overflows might occur in recursion?)
### Stretch Goals
- [ ] Encode Node/Point locations using Morton keys
- [ ] Implement dual tree traversal (DTT) for source/target trees 

### Replace pointers with array indices
Instead of having every Node in the tree point to another struct in memory with
a pointer (64 bits, high search latency because the memory accesses are not contiguous), 
use an unsigned 32-bit integer to manual index an array of Nodes. This allows the array
to be continuous in memory. Additionally, this array should provide the guarrantee that
all 8 child nodes of a given parent are continuous in memory. This will likely help keep 
peer child nodes in lower-level cache more optimally, but may not be helpful for depth-first search.

### Replace recursion with explicit traversal

#### Recursion on a single tree
1. Given a target point `p`, start at the root node of the source tree:   
```c
traverse(root, p);
```
2. Search each of the children of `root`: 
```c
for (uint8_t i=0; i<8; i++) {
    traverse(root.children[i], p);
}
```
3. If the stopping criterion is reached (defined as `theta`), perform calculation; else, continue recursion. Full function:

```c
void traverse(Node node, Point p) {
    if (far_enough_away(node, p, theta)) {
        compute_node_interaction(node, p); 
    }
    else {
        if (has_children(node)) { traverse_children(node, p); }
        else { compute_direct_interaction(node, p); }
    }
}

void traverse_children(Node node, Point p) {
    for (uint8_t i=0; i<8; i++) { 
        traverse(node.children[i], p);
    }
}
```
In this case, the only difference between `compute_node_interaction()` and `compute_direct_interaction()` is that the first function pulls the current density moment from the node (summed from the children) and the second pulls it from a leaf point on the node. This could be further simplified by moving the data from the point into the node and thereby needing only the first function.

#### Explicit Traversal
Explicit traversal means keeping traversing the tree in an orderly fashion or using a stack. 

I don't really want to implement a stack, so we'll do a depth-frst left-right search:
1. Start at the root node  
2. Search the children from left-right
3. If the distance criteria is met, compute interaction; else, descend
4. Once computation of the interaction is complete, ascend one level and look to the next child
5. Ascend from each level once the 8th child is searched; terminate when the pointer is to the root node

```c
typedef struct Target {
    double x, y, z, Bx, By, Bz;
} Point; 

typedef struct Source { 
    double x, y, z, vol, Jx, Jy, Jz;
} Source;

typedef struct Node { 
    double cx, cy, cz, span; 
    int parent, first_child, position;
    Source source; 
} Node; 

typedef struct Octree {
    Node *nodes; 
    Source *sources;
} Octree; 

void traverse(Octree tree, Point p) {

    int node_index = 1;
    Node node = tree.nodes[node_index];
    int parent_index = node.parent;
    int child_offset = 0; 

    while (true) {
        // What node are we looking at? 
        int child_index = node.first_child + child_offset;
        Node child = tree.nodes[child_index];

        // Do we calculate the interaction or descend?
        if (far_enough_away(child, p, theta)) { 
            compute_node_interaction(child, p); 
            child_offset++;
        }
        else {
            if (has_children(child)) { 
                node = child; 
                child_offset = 0; 
            }
            else { 
                compute_direct_interaction(child, p); 
                child_offset++; 
            }
        }

        // Do we ascend? 
        if (child_offset == 8) { 
            // Move up
        }

        if ((node.position == 7)) {
            if (parent_index == 0) { break;} 
            else {
                node_index = parent_index; 
            }
        }
        
    }
}
```

### Dual tree traversal 
The only way I can think of doing DTT is through recursion - I'm not familiar with the stack based approaches to doing this. 

The Biot Savart calculation process using dual tree traversal should look something like the following:
1. Create the source tree 
2. Create the target tree (this may be identical to the source tree for self fields)
3. Partition work amongst threads; for a naive implementation and assuming a balanced octree, this should be done at the first or second level of the tree. up to 8 threads are naturally supported without much trouble 
4. Traverse the two trees using recursion: subdivide the larger of the two nodes, descend

## References to Look At
https://iopscience.iop.org/article/10.1088/1361-6668/aae957/pdf 
https://ieeexplore.ieee.org/document/4202642
https://ieeexplore.ieee.org/document/5754906?arnumber=5754906
https://math.nyu.edu/~greengar/shortcourse_fmm.pdf
https://people.eecs.berkeley.edu/~demmel/cs267/lecture27/lecture27.html
https://people.eecs.berkeley.edu/~demmel/cs267/lecture26/lecture26.html
https://www.cs.princeton.edu/courses/archive/fall04/cos126/assignments/barnes-hut.html
https://beltoforion.de/en/barnes-hut-galaxy-simulator/
https://home.ifa.hawaii.edu/users/barnes/software.html