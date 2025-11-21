#include <stdlib.h> 
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "../include/internal/biotsavart.h"
#include "../include/internal/octree.h"
#include "../include/internal/utils.h"


// Constants for later calculations
const double ONE_OVER_FOUR_THIRDS_PI = 1 / (4.0 * PI / 3.0 );
const double ONE_THIRD = 1.0/3.0;


// --- 
// Data structures
// ---


// Biot Savart source point (from a finite element)
// 6 * 8 bytes = 48 bytes per object; unfortunately this will use an entire cache line
typedef struct Point {
    double x, y, z;             // [m] centroid location vector
    double vJx, vJy, vJz;       // [A-m] current density 'moment', i.e. vol*J
} Point;


// A node, which either can be a leaf or contain eight subnodes 
//
// if a node/point is on the boundary of a subdivision, then it is placed on the 
// positive side, i.e. cx+halfwidth > x >= cx
// Subdivisions are handled in counterclockwise order, + to -, 
// i.e. children[0] is +X,+Y,+Z children[1] is -X,+Y,+Z, children[4] is +X,+Y,-Z, etc.
typedef struct Node {
    double cx, cy, cz;              // centroid in 3D space
    double halfwidth;               // extent; volume of node = (2*halfwidth)^2
    struct Node *children[8];       // null if leaf (no children)
    Point *point;                   // null if node (no point)
    double vJx, vJy, vJz;           // volume*J "moment" 
} Node; 


// --- 
// Tree construction functions
// ---


// Convenience function to create a new source point
// (generally unused)
Point *new_point(double x, double y, double z, double vol, double Jx, double Jy, double Jz) {
    Point *point = malloc(sizeof(Point)); 
    point->x = x; point->y = y; point->z = z;  
    point->vJx = vol*Jx; point->vJy = vol*Jy; point->vJz = vol*Jz; 
    return point;
}


Point *points_from_elements(
    const double *restrict centx, const double *restrict centy, const double *restrict centz, 
    const double *restrict vol, 
    const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n
) {
    Point *points = calloc(n, sizeof(Point));
    for (size_t i=0; i<n; i++) {
        points[i].x = centx[i];
        points[i].y = centy[i];
        points[i].z = centz[i];
        points[i].vJx = vol[i]*Jx[i];
        points[i].vJy = vol[i]*Jy[i];
        points[i].vJz = vol[i]*Jz[i];
    }
    return points;
}


double span_from_coords(
    const double *restrict x, const double *restrict y, const double *restrict z, 
    size_t n
) {
    // returns the side length of the largest enclosing cube for all source pts

    double xmin = min(x, n); double xmax = max(x, n); 
    double ymin = min(y, n); double ymax = max(y, n); 
    double zmin = min(z, n); double zmax = max(z, n); 

    // TODO: make bounding box non-cube
    double xrange = xmax - xmin; 
    double yrange = ymax - ymin; 
    double zrange = zmax - zmin; 

    double span;
    if ((xrange >= yrange) && (xrange >= zrange)) { span = xrange; }
    else if ((yrange >= xrange) && (yrange >= zrange)) { span = yrange; }
    else { span = zrange; }

    return span;
}


// Check if the node is a leaf or a branch
// Leaf -> has no children (FALSE)
// Branch -> has children (TRUE)
bool has_children(Node *node) {
    if (node == NULL) {
        return false;
    }
    for (size_t i=0; i<8; i++) {
        if (node->children[i] != NULL) {
            return true;
        }
    }
    return false;
}


// Make the root node given a centroid and span
// Use the first node in the allocation
Node *make_root(NodeAllocation *nodes, double cx, double cy, double cz, double span) {
    Node *root = &nodes->nodes[0];
    root->cx = cx; root->cy = cy; root->cz = cz; root->halfwidth = span; 
    for (size_t i=0; i<8; i++) {
        root->children[i] = NULL; 
    }
    root->point = NULL; 
    return root;
}


Node *root_from_coords(
    NodeAllocation *nodes, 
    const double *restrict x, const double *restrict y, const double *restrict z, 
    size_t n
) {
    double xmin = min(x, n); double xmax = max(x, n); 
    double ymin = min(y, n); double ymax = max(y, n); 
    double zmin = min(z, n); double zmax = max(z, n); 

    // TODO: make bounding box non-cube (?)
    double xrange = xmax - xmin; 
    double yrange = ymax - ymin; 
    double zrange = zmax - zmin; 
    double cx = (xmax + xmin) / 2.0; 
    double cy = (ymax + ymin) / 2.0; 
    double cz = (zmax + zmin) / 2.0;

    double span;
    if ((xrange >= yrange) && (xrange >= zrange)) { span = xrange; }
    else if ((yrange >= xrange) && (yrange >= zrange)) { span = yrange; }
    else { span = zrange; }

    Node *root = make_root(nodes, cx, cy, cz, span);
    return root;
}


NodeAllocation *allocate_nodes(int n) {
    Node *nodes = calloc(n, sizeof(Node));
    NodeAllocation *node_allocation = calloc(1, sizeof(NodeAllocation));
    node_allocation->nodes = nodes;  
    node_allocation->capacity = n; 
    node_allocation->current_node = 0;
    return node_allocation;
}


int deallocate_nodes(NodeAllocation *allocation) {
    free(allocation->nodes);
    free(allocation);
    return 0;
}


// Increase the size of the node allocation by copying the old allocation and 
// getting a new block of memory
NodeAllocation *reallocate_nodes(NodeAllocation *old_allocation, int new_size) {
    Node *new_nodes = calloc(new_size, sizeof(Node)); 

    for (int i=0; i<old_allocation->capacity; i++) {
        // Copy from old to new (slow, painful)
        new_nodes[i].cx = old_allocation->nodes[i].cx;
        new_nodes[i].cy = old_allocation->nodes[i].cy;
        new_nodes[i].cz = old_allocation->nodes[i].cz;
        new_nodes[i].halfwidth = old_allocation->nodes[i].halfwidth;
        for (int j=0; j<8; j++) {
            new_nodes[i].children[j] = old_allocation->nodes[i].children[j];
        }
        new_nodes[i].point = old_allocation->nodes[i].point; 
        new_nodes[i].vJx = old_allocation->nodes[i].vJx;
        new_nodes[i].vJy = old_allocation->nodes[i].vJy;
        new_nodes[i].vJz = old_allocation->nodes[i].vJz;
    }

    NodeAllocation *new_allocation = calloc(1, sizeof(NodeAllocation));
    new_allocation->nodes = new_nodes;
    new_allocation->capacity=new_size; 
    new_allocation->current_node =old_allocation->current_node;
    int success = deallocate_nodes(old_allocation);
    return new_allocation;
}


// adds a fresh node to the root node at the specified octant
int add_node(NodeAllocation *nodes, Node *root, Octant octant) {
    Node *node;
    nodes->current_node++;
    if (nodes->current_node >= nodes->capacity-1) {
        nodes = reallocate_nodes(nodes, nodes->capacity*2);
    }
    node = &nodes->nodes[nodes->current_node];

    node->halfwidth = root->halfwidth/2.0; 
    double cx, cy, cz;
    double rcx = root->cx; 
    double rcy = root->cy; 
    double rcz = root->cz;
    double rhw2 = root->halfwidth/2;

    switch (octant) {
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

    node->cx = cx; node->cy = cy; node->cz = cz; 
    // child nodes and points already have null pointers from calloc() 
    root->children[octant] = node;

    return 0; 
}


// where is the point in the node region?
Octant find_octant(Point point, Node *node) {
    double x = point.x; double y = point.y; double z = point.z;
    double cx = node->cx; double cy = node->cy; double cz = node->cz;
    // note this does not do a bounds check to see if the point is within the Node region
    if ( (x>=cx) && (y>=cy) && (z>=cz) ) { return NEU; }
    else if ( (x<cx) && (y>=cy) && (z>=cz) ) { return NWU; }
    else if ( (x<cx) && (y<cy) && (z>=cz) ) { return SWU; }
    else if ( (x>=cx) && (y<cy) && (z>=cz) ) { return SEU; }
    else if ( (x>=cx) && (y>=cy) && (z<cz) ) { return NED; }
    else if ( (x<cx) && (y>=cy) && (z<cz) ) { return NWD; }
    else if ( (x<cx) && (y<cy) && (z<cz) ) { return SWD; }
    else if ( (x>=cx) && (y<cy) && (z<cz) ) { return SED; }
    else { return -1; } // error
}


// is the point within the node's domain?
bool check_in_domain(Point point, Node *node) {
    double dx = fabs(node->cx - point.x);
    double dy = fabs(node->cy - point.y);
    double dz = fabs(node->cz - point.z);

    // TODO: this will create an error if there's a node on the boundary (?)
    if ( (dx <= node->halfwidth) && (dy <= node->halfwidth) && (dz <= node->halfwidth) ) {
        return true;
    }
    else {
        return false;
    }
}


int add_points(
    NodeAllocation *nodes, Node *root, Point *points, 
    size_t npts, size_t start_point, size_t end_point
) {

    for (size_t i=start_point; i<end_point; i++) {
        Point *point = &points[i];

        // check that this is a valid point for this domain
        if (!check_in_domain(*point, root) ) {
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
                // both the existing and new point need to go into octants

                // deal with existing point first 
                Octant octant = find_octant(*root->point, root);
                int success = add_node(nodes, root, octant);
                root->children[octant]->point = root->point;
                root->point = NULL;

                // recursively call this function to deal with the new point
                // should only execute once
                add_points(nodes, root, points, npts, i, i+1); 
            }
            else {
                // bin the new point into an octant
                Octant octant = find_octant(*point, root);

                if ( root->children[octant] == NULL ) {
                    // octant node has not already been defined, make one
                    int success = add_node(nodes, root, octant);
                    root->children[octant]->point = point; 
                }
                else {
                    // Octant exists, may need to convert it from a leaf
                    // recursive call again 
                    add_points(nodes, root->children[octant], points, npts, i, i+1);
                }
            }
        }
    }
    
    return 0;
}


// --- 
// Magnetic fields calculations
// --- 


int calculate_moments(Node *root) {

    // first need to recursively sum the moments for each of the children
    if (has_children(root)) {
        for (size_t i=0; i<8; i++) {
            if (root->children[i] != NULL) {
                int success = calculate_moments(root->children[i]); 
                // Now that we've calculated the moments of the child nodes, add them here
                root->vJx += root->children[i]->vJx;
                root->vJy += root->children[i]->vJy;
                root->vJz += root->children[i]->vJz;  
            }
        }
    }
    else {
        // Now we're at the lowest level of the tree
        if ( root->point != NULL) {
            root->vJx += root->point->vJx;
            root->vJy += root->point->vJy;
            root->vJz += root->point->vJz;
        }
    }
    return 0;
}


int print_tree_sum(Node *root) {
    printf("tree sum:  (vJx, vJy, vJz) = (%f, %f, %f)\n", root->vJx, root->vJy, root->vJz);
    return 0;
}


// Compute the contribution of a node to the B-field at a point (x,y,z)
int bfield_node_contribution(
    Node *node, double x, double y, double z, 
    double *Bx, double *By, double *Bz, 
    size_t i, double phi
) {
    double rx = x - node->cx;
    double ry = y - node->cy; 
    double rz = z - node->cz; 
    double rmag = sqrt(rx*rx + ry*ry + rz*rz);

    double sd = (2*node->halfwidth)/rmag;        // theta < s/d (acceptance)
    int success = 0;
    if ((sd > phi) && has_children(node) ) {
        for (size_t j=0; j<8; j++) {
            if (node->children[j] != NULL) {
                success += bfield_node_contribution(node->children[j], x, y, z, Bx, By, Bz, i, phi);
            }
        }
    }
    else {
        // Biot Savart kernel
        double R = node->halfwidth;
        double rmag3 = rmag*rmag*rmag;
        double inv_rmag3;

        // If inside the source radius, correct based on volume of current
        // density enclosed
        if (rmag > R) {inv_rmag3 = 1/rmag3;} 
        else {inv_rmag3 = powf(rmag / R, 3);}

        // Experimental correction determined via numerical integration of a 
        // cube element (not needed)
        // inv_rmag3 *= 1.04*cos(node->halfwidth/rmag);

        // Calculate cross-product (J x r')
        double jxrpx = node->vJy*rz - node->vJz*ry; 
        double jxrpy = node->vJz*rx - node->vJx*rz; 
        double jxrpz = node->vJx*ry - node->vJy*rx;

        // Compute contribution to field 
        Bx[i] += MU04PI * jxrpx * inv_rmag3;
        By[i] += MU04PI * jxrpy * inv_rmag3;
        Bz[i] += MU04PI * jxrpz * inv_rmag3;
    }

    return success;
}