#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double drand(double min, double max) {
    return min + (max-min) * rand() / RAND_MAX;
}

double *vrand(size_t n, double min, double max) {
    double *values = calloc(n, sizeof(double));
    for (size_t i=0; i<n; i++) {
        values[i] = drand(min, max);
    }
    return values; 
}

// Get the maximum value in an array
double max(const double *restrict array, size_t n) {
    double current_max = array[0]; 

    for (size_t i=1; i<n; i++) {
        if (array[i] > current_max) {
            current_max = array[i];
        }
    }
    return current_max;
}

// Get the minimum value in an array
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


// return the root mean square error of `y` measured against `x`
// Note: comparisons where x[i] are zero are ignored
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