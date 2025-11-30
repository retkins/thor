/* Linear dual-tree octree 
*/

#include <stdint.h>         // we use uint32_t throughout because we don't need more than 4B elements in the tree...

// Contains all source points for the source tree
typedef struct Sources {
    uint32_t n;
    double *x, *y, *z; 
    double *vol; 
    double *Jx, *Jy, *Jz;
} Sources; 

// Contains all output points at which the fields will be calculated 
typedef struct Targets {
    uint32_t n; 
    double *x, *y, *z; 
    double *Bx, *By, *Bz
} Targets;

// an index `i` accesses all of these arrays for info on the same node
typedef struct SourceOctree {
    uint32_t n; 
    double size;              // length of the bounding box
    uint32_t *level;            // also defines the size of the grid (grid_size = size/2^level[i])
    double *x, *y, *z;          // centroids of the grid cell (node) in 3D space 
    double *cx, *cy, *cz;       // 'center of mass' of the grid cell (node), weighted by all Sources it contains
    double *vJx, *vJy, *vJz;    // current density-moment of the grid cell (sum of all nodes) 
    uint32_t *parent;           // index of the parent 
    uint32_t *first_child;      // index of the first child (left-most in array/2d view)
    uint32_t *next_sibling;     // index of the sibling to the child's right (MAXINT if none)
    uint32_t *leaf_source;      // if the node is a leaf, links to the Source 
} SourceOctree;

// basically same layout, but vJ is swapped for B
typedef struct TargetOctree {
    uint32_t n;
    double size;                // bounding box size
    uint32_t *level;
    double *x, *y, *z;          // actual coordinates of the grid cell centroid
    double *cx, *cy, *cz;       // centroid of the grid cell as an average of the centroids of all enclosed Source points
    double *Bx, *By, *Bz;       // magnetic field contributions at each level of the grid cell 
    uint32_t *parent; 
    uint32_t *first_child; 
    uint32_t *next_sibling; 
    uint32_t *leaf_target; 
} TargetOctree; 



