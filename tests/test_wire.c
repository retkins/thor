/*  Test that the wire calculation is accurate

    Not currently working...
*/


#include "../src/biotsavart.h"
#include <stdio.h>



int main() {

    const size_t n_sources =2000e4; 
    const size_t n_targets = 1;

    // Solid conductor with centroid at the origin, along the y axis
    double radius = .5; 
    double unit_length = 2*radius;
    double area = unit_length*unit_length;
    double length = n_sources*unit_length; 
    double volume = length*area;
    double unit_volume = volume/n_sources;
    unit_volume = unit_length*unit_length*unit_length;

    double current = 1e3; 
    double Jmag = current/area;
    double y0 = -length/2;
    double *centx = calloc(n_sources, sizeof(double)); 
    double *centy = calloc(n_sources, sizeof(double)); 
    double *centz = calloc(n_sources, sizeof(double));
    double *vol = calloc(n_sources, sizeof(double));
    double *Jx = calloc(n_sources, sizeof(double));
    double *Jy = calloc(n_sources, sizeof(double));
    double *Jz = calloc(n_sources, sizeof(double));
    double dy = length/n_sources;
    double _y = y0;
    _y = 0;
    for (size_t i=0; i<n_sources; i++) {
        centy[i] = _y; 
        _y += dy; 
        Jy[i] = Jmag;
        vol[i] = unit_volume;
    }
    double _x = 2; 
    double *x = calloc(n_targets, sizeof(double)); 
    double *y = calloc(n_targets, sizeof(double)); 
    double *z = calloc(n_targets, sizeof(double)); 
    for (size_t i=0; i<n_targets; i++) {
        x[i] = _x; 
    }
    double *Bx = calloc(n_targets, sizeof(double)); 
    double *By = calloc(n_targets, sizeof(double)); 
    double *Bz = calloc(n_targets, sizeof(double)); 

    int success = bfield_wire(centx, centy, centz, vol, Jx, Jy, Jz, n_sources, x, y, z, Bx, By, Bz, n_targets); 
    double expected_value = MU0*Jmag*area/(2*PI*_x); 
            printf("Calculated value: %.6f\n", Bz[0]);
        printf("Expected value: %.6f\n", expected_value);
    if (fabs(expected_value + Bz[0]) < 1e-8) {
        printf("Test \033[32mPASSED\033[0m.\n");
    }
    else {
        printf("Test \033[31mFAILED\033[0m.\n");

    }
}