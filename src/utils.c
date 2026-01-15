#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>


double drand(double min, double max) {
    return min + (max-min) * rand() / RAND_MAX;
}

double *zeros(uint32_t n) {
    double *array = calloc(n, sizeof(double));
    if (!array) {
        fprintf(stderr, "Error in `zeros()`: memory allocation of %i elements failed.\n", n);
        return NULL;
    }
    else {
        return array;
    }
}


double *vrand(size_t n, double min, double max) {
    double *values = calloc(n, sizeof(double));
    for (size_t i=0; i<n; i++) {
        values[i] = drand(min, max);
    }
    return values; 
}


double max(const double *restrict array, size_t n) {
    double current_max = array[0]; 

    for (size_t i=1; i<n; i++) {
        if (array[i] > current_max) {
            current_max = array[i];
        }
    }
    return current_max;
}


double min(const double *restrict array, size_t n) {
    double current_min = array[0]; 

    for (size_t i=1; i<n; i++) {
        if (array[i] < current_min) {
            current_min = array[i];
        }
    }
    return current_min;
}


double test(double x) {
    printf("Hello from C!\n");
    return 2*x;
}


double rms_error(double *restrict y, double *restrict x, size_t n) {

    double error = 0.0; 
    for (size_t i=0; i<n; i++) {
        if (x[i] > 1e-8) {
            double _error = (y[i] - x[i]) / x[i];
            error += _error*_error;
        }
    }

    error /= n; 
    return sqrt(error);
}

static inline void cross(double a1, double a2, double a3, double b1, double b2, double b3, double *c1, double *c2, double *c3) {
    *c1 = a2*b3 - a3*b2; 
    *c2 = a3*b1 - a1*b3; 
    *c3 = a1*b2 - a2*b1; 
}

static inline double norm(double a1, double a2, double a3) {
    return fabs(a1*a1 + a2*a2 + a3*a3);
}