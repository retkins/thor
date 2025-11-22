/*  Setup dicrete problems useful for running tests or evaluating standard 
    electromagnetic conditions.

    Define the following:
    - Current carrying circular loop 
    - Simple solenoid
    - Dense solenoid
*/

#ifndef PROBLEMS_H
#define PROBLEMS_H 

#include <stdint.h>

// ---
// Simple current-carrying circular loop 
// ---

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

// ---
// Simple solenoid
// ---

// Simple solenoid with one element in radial direction
typedef struct Solenoid {
    double radius;              // [m]
    double length;              // [m]
    double current;             // [A]
    uint32_t n;                 // number of points
    double *x, *y, *z;          // [m] location of source points in loop 
    double *vol;                // [m^3] volume of each source point 
    double *Jx, *Jy, *Jz;       // [A/m^2] current density vector at each source point

} Solenoid;

// --- 
// Output arrays
// --- 

// Defines which direction a series of points should be drawn along
typedef enum Direction {X, Y, Z} Direction; 

// A line in 3D space drawn along the specified axis
typedef struct Line {
    uint32_t n;                 // number of points
    double *x, *y, *z;          // [m]
    double *Bx, *By, *Bz;       // [T]
} Line;

// Create a new line along the specified global axis
Line *new_line(Direction dir, double start, double end, uint32_t n);

// Free allocated memory for a Line
void free_line(Line *line);

#endif