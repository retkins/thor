/*  Routines for handling Morton (Z-order) codes
*/

#ifndef MORTON_H
#define MORTON_H

#include <stdint.h>


// Define the bounding box for the space that an octree represents
typedef struct BBox {
    double xc, yc, zc; 
    double side_length;             // TODO optimize by storing 1/side_length
    double xmin, ymin, zmin;        // do we need to store both centroid and min?
} BBox;

// Constructore for BBox
BBox new_bbox(double xc, double yc, double zc, double side_length);

// Encode a series of points into Morton representation
// L is the number of levels in the tree (L<=21 is standard)
uint64_t *encode_points(double *x, double *y, double *z, uint32_t n, uint32_t L, BBox bbox); 

#endif