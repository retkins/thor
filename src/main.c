#include <stdlib.h> 
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

// Define structures to represent the octree

// an end point (leaf), which is a biot savart source element
typedef struct {
    double x, y, z; 
    double vol; 
    double Jx, Jy, Jz;
} Point;

// Convenience function to create a new source point
Point *new_point(double x, double y, double z, double vol, double Jx, double Jy, double Jz) {
    Point *point = malloc(sizeof(Point)); 
    point->x = x; point->y = y; point->z = z; 
    point->vol = vol; 
    point->Jx = Jx; point->Jy = Jy; point->Jz = Jz; 
    return point;
}

// north/south/east/west/up/down
typedef enum Octrant { NEU, NWU, SWU, SEU, NED, NWD, SWD, SED } Octrant; 

// a node, which either can be a leaf or contain eight subnodes 
// if a node/point is on the boundary of a subdivision, then it is placed on the 
// positive side, i.e. cx+halfwidth > x >= cx
// Subdivisions are handled in counterclockwise order, + to -, 
// i.e. children[0] is +X,+Y,+Z children[1] is -X,+Y,+Z, children[4] is +X,+Y,-Z, etc.
typedef struct Node {
    double cx, cy, cz;              // centroid in 3D space
    double halfwidth;               // extent; volume of node = (2*halfwidth)^2
    struct Node *children[8];       // null if leaf (no children)
    Point *point;                   // null if node (no point)
} Node; 

// check if the node is a leaf 
bool has_children(Node *node) {
    for (size_t i=0; i<8; i++) {
        if (node->children[i] != NULL) {
            return true;
        }
    }
    return false;
}

// adds a fresh node to the root node at the specified octrant
int add_node(Node *root, Octrant octrant) {
    Node *node = calloc(1, sizeof(Node)); 
    node->halfwidth = root->halfwidth/2.0; 
    double cx, cy, cz;
    double rcx = root->cx; 
    double rcy = root->cy; 
    double rcz = root->cz;
    double rhw2 = root->halfwidth/2;

    switch (octrant) {
        case NEU: 
            cx = rcx + rhw2; cy = rcy + rhw2; cz = rcz + rhw2; break;
        case NWU: 
            cx = rcx - rhw2; cy = rcy + rhw2; cz = rcz + rhw2; break;
        case SWU: 
            cx = rcx - rhw2; cy = rcy - rhw2; cz = rcz + rhw2; break;
        case SEU: 
            cx = rcx + rhw2; cy = rcy - rhw2; cz = rcz + rhw2; break;
        case NED: 
            cx = rcx + rhw2; cy = rcy + rhw2; cz = rcz - rhw2; break;
        case NWD: 
            cx = rcx - rhw2; cy = rcy + rhw2; cz = rcz - rhw2; break;
        case SWD: 
            cx = rcx - rhw2; cy = rcy - rhw2; cz = rcz - rhw2; break;
        case SED: 
            cx = rcx + rhw2; cy = rcy - rhw2; cz = rcz - rhw2; break;
        default:
            printf("error in add_node()!\n"); break;
    }

    node->cx = cx; node->cy = cy; node->cz = cz; node->point=NULL; 
    root->children[octrant] = node;
    return 0; 
}

// where is the point in the node region?
Octrant find_octrant(Point *point, Node *node) {
    double x = point->x; double y = point->y; double z = point->z;
    // note this does not do a bounds check to see if the point is within the Node region
    if ( (x>=node->cx) && (y>=node->cy) && (z>=node->cz) ) { return NEU; }
    else if ( (x<node->cx) && (y>=node->cy) && (z>=node->cz) ) { return NWU; }
    else if ( (x<node->cx) && (y<node->cy) && (z>=node->cz) ) { return SWU; }
    else if ( (x>=node->cx) && (y<node->cy) && (z>=node->cz) ) { return SEU; }
    else if ( (x>=node->cx) && (y>=node->cy) && (z<node->cz) ) { return NED; }
    else if ( (x<node->cx) && (y>=node->cy) && (z<node->cz) ) { return NWD; }
    else if ( (x<node->cx) && (y<node->cy) && (z<node->cz) ) { return SWD; }
    else if ( (x>=node->cx) && (y<node->cy) && (z<node->cz) ) { return SED; }
    else { return -1; } // error
}

// is the point within the node's domain?
bool check_in_domain(Point *point, Node *node) {
    double dx = fabs(node->cx - point->x);
    double dy = fabs(node->cy - point->y);
    double dz = fabs(node->cz - point->z);

    // TODO: this will create an error if there's a node on the boundary
    if ( (dx <= node->halfwidth) && (dy <= node->halfwidth) && (dz <= node->halfwidth) ) {
        return true;
    }
    else {
        return false;
    }
}

// add Points to the tree
int add_points(Node *root, Point **points, size_t npts, size_t start_point, size_t end_point) {

    // remember this is a recursive function...
    for (size_t i=start_point; i<end_point; i++) {
        Point *point = points[i];

        // check that this is a valid point for this domain
        if (!check_in_domain(point, root) ) {
            // It's outside the boundary, we skip it
            // note that due to recursion we may do many domain checks here; 
            //  room to optimize
            continue;
        }
        else {
            // its inside the boundary, so we need to find a home for it

            if (!has_children(root) && root->point == NULL) {
                // the root has no children and no point; make it a leaf
                root->point = point;
            }
            else if (!has_children(root) && root->point != NULL) {
                // there are no children but there is an existing point
                // both the existing and new point need to go into octrants

                // deal with existing point first 
                Octrant octrant = find_octrant(root->point, root);
                int success = add_node(root, octrant);
                root->children[octrant]->point = root->point;
                root->point = NULL;

                // recursively call this function to deal with the new point
                add_points(root, points, npts, i, i+1); // should only execute once
            }
            else {
                // bin the new point into an octrant
                Octrant octrant = find_octrant(point, root);

                if ( root->children[octrant] == NULL ) {
                    // octrant node has not already been defined, make one
                    int success = add_node(root, octrant);
                    root->children[octrant]->point = point;
                }
                else {
                    // Octrant exists, may need to convert it from a leaf
                    // recursive call again 
                    add_points(root->children[octrant], points, npts, i, i+1);
                    // printf("test\n");
                }
            }
        }
    }
    
    return 0;
}

// test data 
// const size_t npts = 3;
// const double coords[3*npts] = { 
//     0.5, 0.5, 0.0, 
//     0.4, 0.5, 0.0, 
//     -0.6, -0.2, 0.0
// };

double drand(double min, double max) {
    return min + (max-min) * rand() / RAND_MAX;
}

double *vrand(size_t n) {
    double *values = calloc(n, sizeof(double));
    for (size_t i=0; i<n; i++) {
        values[i] = drand(-1.0, 1.0);
    }
    return values; 
}

// for testing
Point **points_from_array(double* coords, size_t npts) {
    Point **points = malloc(npts*sizeof(Point*));
    size_t ncoords = (size_t)3*npts;
    size_t j=0;
    for (size_t i=0; i<npts; i++) {
        points[i] = malloc(sizeof(Point));
        points[i]->x = coords[j++]; 
        points[i]->y = coords[j++];
        points[i]->z = coords[j++];
    }

    return points;
}


// TODO: none of the memory is freed, relying on the OS once the program exits..
int main() {

    // For testing, assume a root node at origin with half width 1m
    Node *root = malloc(sizeof(Node));
    root->cx = 0.0; root->cy = 0.0; root->cz = 0.0; 
    root->halfwidth = 1.0; 

    const size_t npts = 1000000; 
    const size_t ncoords = 3*npts;
    double *coords = vrand(ncoords); 

    // test data
    Point **points = points_from_array(coords, npts); 

    // test to see if it runs...
    clock_t start_time = clock(); 
    int success = add_points(root, points, npts, 0, npts);
    clock_t end_time = clock(); 

    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("npts: %li, elapsed: %f s\n", npts, elapsed_time);

    if (success == 0) {
        printf("Program completely successfully.\n");
    }
    else {
        printf("Unknown error.\n");
    }

    return 0;
}

// <0.5s to build the octree for 1M points on an M1 Pro (2021)... not bad?