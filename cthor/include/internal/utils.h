

#ifndef UTILS_H 
#define UTILS_H
#include <stdlib.h>
#include <stdint.h>

// Create a random double value in the range (min, max)
double drand(double min, double max);

// Create an array of doubles on the heap (zeros) 
double *zeros(uint32_t n); 

// Create a random double array on the heap
double *vrand(size_t n, double min, double max);

// Get the minimum value in an array
double min(const double *restrict array, size_t n); 

// Get the maximum value in an array
double max(const double *restrict array, size_t n);

// Test for calling from julia
double test(double x);

// Return the root mean square error of `y` measured against `x`
// Note: comparisons where x[i] are zero are ignored
// 
// RMSE = sqrt(1/N * sum_i( (yi - xi)/xi)^2 )
double rms_error(double *restrict y, double *restrict x, size_t n);

// In-place cross-product 
// a x b = c
static inline void cross(double a1, double a2, double a3, double b1, double b2, double b3, double *c1, double *c2, double *c3);
static inline double norm(double a1, double a2, double a3);
#endif