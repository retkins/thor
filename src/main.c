#include <stdlib.h> 
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "utils.h"
#include "octree.h"
#include "biotsavart.h"


// TODO: none of the memory is freed, relying on the OS once the program exits..
int main_old() {

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

    if (success == 0) {
        printf("Program completely successfully.\n");
    }
    else {
        printf("Unknown error.\n");
    }

    return 0;
}

// <0.5s to build the octree for 1M points on an M1 Pro (2021)... not bad?


int main() {
    // Magnetic flux density calculations

    size_t n = 10000; 
    size_t m = 10000; 
    double phi = 0.1;
    printf("Sources x obs points: (%li, %li)\n", n, m);
    printf("phi = %f\n", phi);

    double *x = calloc(m, sizeof(double));
    double *y = calloc(m, sizeof(double));
    double *z = calloc(m, sizeof(double));
    double *Bx = calloc(m, sizeof(double));
    double *By = calloc(m, sizeof(double));
    double *Bz = calloc(m, sizeof(double));
    double *_Bx = calloc(m, sizeof(double));
    double *_By = calloc(m, sizeof(double));
    double *_Bz = calloc(m, sizeof(double));
    for (size_t i=0; i<m; i++) {
        // magic numbers for now 
        x[i] = 1.0; // drand(-1,1); 
        y[i] = 0.0; // drand(-1,1); 
        z[i] = 0.25; // drand(-1,1);
    }

    double *vol = calloc(n, sizeof(double));
    double *centx = calloc(n, sizeof(double));
    double *centy = calloc(n, sizeof(double));
    double *centz = calloc(n, sizeof(double));
    double *Jx = calloc(n, sizeof(double));
    double *Jy = calloc(n, sizeof(double));
    double *Jz = calloc(n, sizeof(double));

    double theta = 0; 
    double dtheta = 2*PI/n;
    double radius = 1.0;
    double side_length = 1e-2;
    double area = side_length*side_length;
    double volume = area*2*PI*radius/n;
    double current = 1e6;
    double current_density = current/area;
    for (size_t i=0; i<n; i++) {
        // magic numbers for now 
        vol[i] = volume; 
        centx[i] = radius*cos(theta);
        centy[i] = radius*sin(theta); 
        centz[i] = 0.0;
        Jx[i] = -current_density*sin(theta);
        Jy[i] = current_density*cos(theta);
        Jz[i] = 0.0;
        theta += dtheta;
    }

    clock_t start_time = clock();
    int success = bfield_direct(centx, centy, centz, vol, Jx, Jy, Jz, n, x, y, z, m, Bx, By, Bz, 1);
    clock_t end_time = clock(); 
    double naive_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Elapsed time for naive (N^2):       %f s\n", naive_time);
    start_time = clock();
    success = bfield_octree(centx, centy, centz, vol, Jx, Jy, Jz, n, x, y, z, m, _Bx, _By, _Bz, 1, phi);
    end_time = clock(); 
    double octree_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Elapsed time for octree (N log(N)): %f s\n", octree_time);
    printf("Speedup = %.3fx\n", naive_time/octree_time);
    printf("first point: \n"); 
    printf("naive:  (Bx, By, Bz) = (%f, %f, %f)\n", Bx[0], By[0], Bz[0]);
    printf("octree: (Bx, By, Bz) = (%f, %f, %f)\n", _Bx[0], _By[0], _Bz[0]);
    double dBx = 0.0;
    double dBy = 0.0;
    double dBz = 0.0;
    for (size_t i=0; i<m; i++) { 
        dBx += pow(Bx[i] - _Bx[i],2); 
        dBy += pow(By[i] - _By[i], 2); 
        dBz += pow(Bz[i] - _Bz[i],2);
    }
    dBx = sqrt(dBx/(double)n); 
    dBy = sqrt(dBy/(double)n); 
    dBz = sqrt(dBz/(double)n); 

    
    printf("x error = %.3f %%\n", dBx*100);
    printf("y error = %.3f %%\n", dBy*100);
    printf("z error = %.3f %%\n", dBz*100);

    return 0;

}