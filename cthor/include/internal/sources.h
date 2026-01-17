/*  Setup dicrete problems useful for running tests or evaluating standard 
    electromagnetic conditions.

    Define the following:
    - Current carrying circular loop 
    - Simple solenoid
    - Dense solenoid
*/

#ifndef SOURCES_H
#define SOURCES_H 

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

// Simple solenoid with one element in radial direction, centered about the z-axis
typedef struct Solenoid {
    double radius;                  // [m]
    double length;                  // [m]
    double zcentroid;               // [m]
    double current;                 // [A]
    double turns_per_unit_length;   // [turns/m]
    uint32_t n;                 // number of points
    double *x, *y, *z;          // [m] location of source points in loop 
    double *vol;                // [m^3] volume of each source point 
    double *Jx, *Jy, *Jz;       // [A/m^2] current density vector at each source point

} Solenoid;

// Allocate memory for a new Solenoid
Solenoid *new_solenoid(double radius, double length, double zcentroid, 
                        double current, double turns_per_unit_length, uint32_t n);

// Free memory for a Solenoid
void free_solenoid(Solenoid *solenoid);

// --- 
// Dense Solenoid
// ---

// Dense solenoid with multiple elements in radial direction, 
// centered about the z-axis
typedef struct DenseSolenoid {
    double inner_radius;            // [m]
    double outer_radius;            // [m]
    double length;                  // [m]
    double zcentroid;               // [m]
    double current;                 // [A]
    double turns_per_unit_length;   // [turns/m]
    double element_size;            // [m]
    uint32_t n;                     // number of points
    double *x, *y, *z;          // [m] location of source points in loop 
    double *vol;                // [m^3] volume of each source point 
    double *Jx, *Jy, *Jz;       // [A/m^2] current density vector at each source point

} DenseSolenoid;

// Allocate memory for a new DenseSolenoid
DenseSolenoid *new_dense_solenoid(
    double inner_radius, double outer_radius, double length, double zcentroid, 
    double current, double element_size
);

// Deallocate memory for a DenseSolenoid
void free_dense_solenoid(DenseSolenoid *dense_solenoid);


#endif