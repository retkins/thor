#include <inttypes.h>

#include "refactor_octree.h"



typedef struct BasicLeaf {
    uint32_t parent; 
    uint32_t source; 
    double vJx, vJy, vJz; 
} BasicLeaf; 

typedef struct BasicBranch {
    uint32_t parent; 
    uint32_t first_child; 
    double vJx, vJy, vJz;
} BasicBranch;

typedef struct Source {
    double x, y, z;         // location in 3D space 
    double vol;             // volume of the source 
    double Jx, Jy, Jz;      // Current density vector 
} Source; 

typedef struct BasicNode {
    uint32_t parent;        // Parent node in the tree 
    uint32_t first_child;   // Index of first child in the array of nodes 
    double cx, cy, cz;      // Centroid of the space that the Node represents
    double halfspan;        // half-width of the cube of space the Node represents
    double vJx, vJy, vJz;   // Current-density moment
    uint32_t source;        // Index of source if Node is a Leaf; else = UINT32_MAX
} BasicNode; 

typedef struct Octree {
    uint32_t n_sources; 
    uint32_t n_nodes;
    uint32_t max_n_nodes;

    Source *sources; 
    BasicNode *nodes;
} Octree;