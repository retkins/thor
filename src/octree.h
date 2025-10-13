#include <stdlib.h> 
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#ifndef OCTREE_H 
#define OCTREE_H 

// an end point (leaf), which is a biot savart source element
typedef struct Point Point;

// Convenience function to create a new source point
Point *new_point(double x, double y, double z, double vol, double Jx, double Jy, double Jz);

double span_from_coords(const double *restrict x, const double *restrict y, const double *restrict z, size_t n);

// north/south/east/west/up/down
typedef enum Octrant { NEU, NWU, SWU, SEU, NED, NWD, SWD, SED } Octrant; 

typedef struct Node Node;

// add Points to the tree
int add_points(Node *root, Point **points, size_t npts, size_t start_point, size_t end_point);

//// calculate the current density-moment for each node in the tree
int calculate_moments(Node *root);

Node *make_root(double cx, double cy, double cz, double span);
Node *root_from_coords(const double *restrict x, const double *restrict y, const double *restrict z, size_t n);
int print_tree_sum(Node *root);
void print_points_sum(Point **points, size_t npts);

Point **points_from_array(double* coords, size_t npts);
#endif 