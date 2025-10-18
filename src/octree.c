#include <stdlib.h> 
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "octree.h"
#include "utils.h"
#include "biotsavart.h"

// Define structures to represent the octree

// an end point (leaf), which is a biot savart source element
typedef struct Point {
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

Point **points_from_elements(
    const double *restrict centx, const double *restrict centy, const double *restrict centz, 
    const double *restrict vol, 
    const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n
) {
    Point **points = calloc(n, sizeof(Point*));
    for (size_t i=0; i<n; i++) {
        points[i] = calloc(1, sizeof(Point));
        points[i]->x = centx[i];
        points[i]->y = centy[i];
        points[i]->z = centz[i];
        points[i]->vol = vol[i]; 
        points[i]->Jx = Jx[i];
        points[i]->Jy = Jy[i];
        points[i]->Jz = Jz[i];
    }
    return points;
}

double span_from_coords(const double *restrict x, const double *restrict y, const double *restrict z, size_t n) {
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

// north/south/east/west/up/down
//typedef enum Octrant { NEU, NWU, SWU, SEU, NED, NWD, SWD, SED } Octrant; 

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
    double vJx, vJy, vJz;           // volume*J "moment" 
} Node; 

// check if the node is a leaf 
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

Node *make_root(double cx, double cy, double cz, double span) {
    Node *root = calloc(1, sizeof(Node));
    root->cx = cx; root->cy = cy; root->cz = cz; root->halfwidth = span; 
    for (size_t i=0; i<8; i++) {
        root->children[i] = NULL; 
    }
    root->point = NULL; 
    return root;
}

Node *root_from_coords(const double *restrict x, const double *restrict y, const double *restrict z, size_t n) {
    double xmin = min(x, n); double xmax = max(x, n); 
    double ymin = min(y, n); double ymax = max(y, n); 
    double zmin = min(z, n); double zmax = max(z, n); 

    // TODO: make bounding box non-cube
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

    Node *root = make_root(cx, cy, cz, span);

    return root;
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
    for (size_t i=0; i<8; i++) {
        node->children[i] = NULL;
    }
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
                }
            }
        }
    }
    
    return 0;
}

// calculate the current density-moment for each node in the tree
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
            root->vJx += root->point->vol*root->point->Jx;
            root->vJy += root->point->vol*root->point->Jy;
            root->vJz += root->point->vol*root->point->Jz;
        }
    }
    return 0;
}

int print_tree_sum(Node *root) {
    printf("tree sum:  (vJx, vJy, vJz) = (%f, %f, %f)\n", root->vJx, root->vJy, root->vJz);
    return 0;
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
        points[i]->vol = 1e-4;  // drand(1e-8, 1e-4);
        points[i]->Jx = 1e6;    // drand(-1e8,1e8);
        points[i]->Jy = 1e6;    // drand(-1e8,1e8);
        points[i]->Jz = 1e6;    // drand(-1e8,1e8);
    }

    return points;
}

void print_points_sum(Point **points, size_t npts) {
    double vJx = 0; double vJy = 0; double vJz = 0; 
    for (size_t i=0; i<npts; i++) {
        vJx += points[i]->vol*points[i]->Jx;
        vJy += points[i]->vol*points[i]->Jy;
        vJz += points[i]->vol*points[i]->Jz;
    }
    printf("array sum: (vJx, vJy, vJz) = (%f, %f, %f)\n", vJx, vJy, vJz);
}


int bfield_node_contribution(Node *node, double x, double y, double z, double *Bx, double *By, double *Bz, size_t i, double phi) {
    double cx = x - node->cx;
    double cy = y - node->cy; 
    double cz = z - node->cz; 
    double rmag = sqrt(cx*cx + cy*cy + cz*cz);

    double sd = (2*node->halfwidth)/rmag;        // s/d
    int success = 0;
    if ((sd > phi) && has_children(node) ) {
        for (size_t j=0; j<8; j++) {
            if (node->children[j] != NULL) {
                success += bfield_node_contribution(node->children[j], x, y, z, Bx, By, Bz, i, phi);
            }
        }
    }
    else {
        double centx, centy, centz;
        if (node->point == NULL) {
            centx = node->cx; centy = node->cy; centz = node->cz;
        }
        else {
            centx = node->point->x; centy = node->point->y; centz = node->point->z;
        }

        // Biot Savart kernel
        double rx = x - centx;
        double ry = y - centy; 
        double rz = z - centz; 
        rmag = rx*rx + ry*ry + rz*rz;
        rmag = sqrt(rmag); 
        double rmag3 = rmag*rmag*rmag;
        double inv_rmag3 = 1/rmag3;
        // inv_rmag3 *= 1.04*cos(node->halfwidth/rmag);

        // Calculate cross-product 
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