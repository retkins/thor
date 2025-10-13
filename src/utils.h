#include <stdlib.h>

#ifndef UTILS_H 
#define UTILS_H

// create a random double value in the range (min, max)
double drand(double min, double max);

// create a random double array on the heap
double *vrand(size_t n, double min, double max);

// Get the minimum value in an array
double min(const double *restrict array, size_t n); 

// Get the maximum value in an array
double max(const double *restrict array, size_t n);

#endif