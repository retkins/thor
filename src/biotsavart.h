#include <math.h>


#ifndef BIOTSAVART_H 
#define BIOTSAVART_H

// Constants for readability 
static const double PI = 3.141592654; 
static const double MU0 = 4*PI*1e-7; 
static const double MU04PI = 1e-7;

// Naive integration of the Biot Savart law for point sources on arbitrary obs pts
// B = mu0/4pi * vol * J x r' / |r'|^3
int bfield_naive(
    const double *restrict centx, const double *restrict centy, const double *restrict centz, 
    const double *restrict vol, const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n, // sources
    const double *restrict x, const double *restrict y, const double *restrict z, 
    size_t m, // obs pts
    double *restrict Bx, double *restrict By, double *restrict Bz, 
    int nthreads
);

// Naive integration of the Biot Savart law for point sources and self-fields
// B = mu0/4pi * vol * J x r' / |r'|^3
int bfield_self_naive(
    const double *restrict centx, const double *restrict centy, const double *restrict centz, 
    const double *restrict vol, const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n,
    double *restrict Bx, double *restrict By, double *restrict Bz
);

// Naive integration of the Biot Savart law for point sources and self_fields using octree method
// B = mu0/4pi * vol * J x r' / |r'|^3
int bfield_self_naive_octree(
    const double *restrict centx, const double *restrict centy, const double *restrict centz, 
    const double *restrict vol, const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n, // sources
    const double *restrict x, const double *restrict y, const double *restrict z, 
    size_t m, // obs pts
    double *restrict Bx, double *restrict By, double *restrict Bz, 
    int nthreads, 
    double phi
);

#endif