#ifndef PROBLEMS_H
#define PROBLEMS_H 

#include <stdint.h>

// Define a current-carrying loop in 3d space to be used as a biot-savart source
// Current loop is aligned along the z-axis and is in the xy plane
typedef struct Loop {
    double radius;              // [m]
    double current;             // [A]
    uint32_t n;                 // number of points
    double *x, *y, *z;          // [m] location of source points in loop 
    double *vol;                // [m^3] volume of each source point 
    double *Jx, *Jy, *Jz;       // [A/m^2] current density vector at each source point
} Loop;

// Create a new Loop (on the heap) 
Loop *new_loop(double radius, double current, uint32_t n);

// Free memory allocated for a Loop
void free_loop(Loop *loop);

#endif