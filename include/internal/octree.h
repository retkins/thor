#include <stdlib.h> 
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#ifndef OCTREE_H 
#define OCTREE_H 

// An end point (leaf), which is a Biot Savart source element
typedef struct Point Point;

// A node in the octree
// Represents a recursive series of smaller nodes within it
typedef struct Node Node;

// Define the locations of each octant in a Node
// north/south/east/west/up/down
typedef enum Octant { NEU, NWU, SWU, SEU, NED, NWD, SWD, SED } Octant; 

// "Arena-style" allocation for Nodes (really, just a fancy array)
// TODO: make this opaque, for some reason biotsavart.c doesn't like it when the 
// fields are defined within the .c file
typedef struct NodeAllocation {
    Node *nodes; 
    int capacity;
    int current_node;
} NodeAllocation;

// Determine the span of the root node from the source point coordinates
// the span is the total size (i.e. side length) of the root node, a cube
// TODO: should this be private inside the .c file?
double span_from_coords(
    const double *restrict x, const double *restrict y, const double *restrict z,
    size_t n
);

// Allocate an array of Point objects that contain the source element data
Point *points_from_elements(
    const double *restrict centx, const double *restrict centy, const double *restrict centz, 
    const double *restrict vol, 
    const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n
);

// Allocate memory for `n` nodes
NodeAllocation *allocate_nodes(int n);

// Free the octree
int deallocate_nodes(NodeAllocation *allocation);

// Add source points to the tree by recursively performing a top-down traversal
int add_points(
    NodeAllocation *nodes, Node *root, Point *points, 
    size_t npts, size_t start_point, size_t end_point
);

// Calculate the current density-moment for each node in the tree (upwards pass)
// The quantity of each child node is added to the parent
// Current-density moment is defined as the vector quantity (vol*J)
int calculate_moments(Node *root);

// Make the root node in the tree from the coordinates of all of the source points
Node *root_from_coords(
    NodeAllocation *nodes, 
    const double *restrict x, const double *restrict y, const double *restrict z, 
    size_t n
);

// Print the sum of the current density moments for the entire tree
// (generally for debug purposes only)
int print_tree_sum(Node *root);

// Compute the magnetic flux density contribution of a node at a target point 
// (x,y,z), which is stored at location `i` in the output B arrays
// `phi` is the angle-opening criteria, which should be <= 0.1 for most applications
// 
// This function is meant to be called recursively, starting at the root node of the tree
void bfield_node_contribution(
    Node *node, double x, double y, double z, 
    double *Bx, double *By, double *Bz, 
    size_t i, double phi
);

#endif 