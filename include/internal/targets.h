/* Define output collections of target points at which fields should be
    calculated
*/

#ifndef TARGETS_H
#define TARGETS_H

#include <stdint.h>

// Defines which cartesian axis a series of points should be drawn along
typedef enum Axis {X, Y, Z} Axis; 

// A line in 3D space drawn along the specified axis
typedef struct Line {
    uint32_t n;                 // number of points
    double *x, *y, *z;          // [m]
    double *Bx, *By, *Bz;       // [T]
} Line;

// Create a new line along the specified global axis
Line *new_line(Axis axis, double start, double end, uint32_t n);

// Free allocated memory for a Line
void free_line(Line *line);

#endif