#include <stdlib.h> 
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "utils.h"
#include "octree.h"
#include "biotsavart.h"


// TODO: none of the memory is freed, relying on the OS once the program exits..
int main() {

    // For testing, assume a root node at origin with half width 1m
    Node *root = make_root(0, 0, 0, 1.0);

    const size_t npts = 1000000; 
    const size_t ncoords = 3*npts;
    double *coords = vrand(ncoords, -1, 1); 

    printf("Running thor with %li npts\n", npts);

    // test data
    Point **points = points_from_array(coords, npts); 

    // test to see if it runs...
    clock_t start_time = clock(); 
    int success = add_points(root, points, npts, 0, npts);
    clock_t end_time = clock(); 

    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Time for tree creation: %f s\n", elapsed_time);

    start_time = clock();
    success = calculate_moments(root);
    end_time = clock();
    elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Time for tree sum:      %f s\n", elapsed_time);
    print_points_sum(points, npts);
    print_tree_sum(root);

    // Magnetic flux density calculations
    size_t m = 1;
    double *x = calloc(m, sizeof(double));
    double *y = calloc(m, sizeof(double));
    double *z = calloc(m, sizeof(double));
    double *Bx = calloc(m, sizeof(double));
    double *By = calloc(m, sizeof(double));
    double *Bz = calloc(m, sizeof(double));
    for (size_t i=0; i<m; i++) {
        // magic numbers for now 
        x[i] = 1.0; 
        y[i] = 1.0; 
        z[i] = 1.0;
    }

    double *vol = calloc(npts, sizeof(double));
    double *centx = calloc(npts, sizeof(double));
    double *centy = calloc(npts, sizeof(double));
    double *centz = calloc(npts, sizeof(double));
    double *Jx = calloc(npts, sizeof(double));
    double *Jy = calloc(npts, sizeof(double));
    double *Jz = calloc(npts, sizeof(double));
    for (size_t i=0; i<npts; i++) {
        // magic numbers for now 
        vol[i] = 1e-4; 
        centx[i] = points[i]
        centy[i] = 1e6; 
        centz[i] = 1e6; 
    }

    success = bfield_naive()

    if (success == 0) {
        printf("Program completely successfully.\n");
    }
    else {
        printf("Unknown error.\n");
    }

    return 0;
}

// <0.5s to build the octree for 1M points on an M1 Pro (2021)... not bad?