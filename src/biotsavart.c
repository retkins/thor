#include <stdlib.h>
#include <time.h>
// #include <omp.h>

#include "biotsavart.h"
#include "octree.h"
#include "utils.h"

// Naive integration of the Biot Savart law for point sources 
// B = mu0/4pi * vol * J x r' / |r'|^3
int bfield_naive(
    const double *restrict centx, const double *restrict centy, const double *restrict centz, 
    const double *restrict vol, const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n, // sources
    const double *restrict x, const double *restrict y, const double *restrict z, 
    size_t m, // obs pts
    double *restrict Bx, double *restrict By, double *restrict Bz, 
    int nthreads
) {

    // Outer loop over sources
    // #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1) reduction(+:Bx[:m], By[:m], Bz[:m])
    for (size_t i=0; i<n; i++) {

        // Inner loop over obs pts 
        for (size_t j=0; j<m; j++) {
            // Calculate r'
            double rx = x[j] - centx[i];
            double ry = y[j] - centy[i]; 
            double rz = z[j] - centz[i]; 
            double rmag = rx*rx + ry*ry + rz*rz;
            rmag = sqrt(rmag); 
            rmag = rmag*rmag*rmag;
            rmag = 1/rmag;

            // Calculate cross-product 
            double jxrpx = Jy[i]*rz - Jz[i]*ry; 
            double jxrpy = Jz[i]*rx - Jx[i]*rz; 
            double jxrpz = Jx[i]*ry - Jy[i]*rx;

            // Compute contribution to field 
            Bx[j] += MU04PI * vol[i] * jxrpx * rmag;
            By[j] += MU04PI * vol[i] * jxrpy * rmag;
            Bz[j] += MU04PI * vol[i] * jxrpz * rmag;
        }
    }

    return 0;
}

// Naive integration of the Biot Savart law for point sources and self_fields
// B = mu0/4pi * vol * J x r' / |r'|^3
int bfield_self_naive(
    const double *restrict centx, const double *restrict centy, const double *restrict centz, 
    const double *restrict vol, const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n, // sources and obs points
    // const double *restrict x, const double *restrict y, const double *restrict z, 
    // size_t m, // obs pts
    double *restrict Bx, double *restrict By, double *restrict Bz
    // int nthreads
) {

    // Outer loop over sources
    // #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1) reduction(+:Bx[:m], By[:m], Bz[:m])
    for (size_t i=0; i<n; i++) {

        // Inner loop over obs pts 
        for (size_t j=0; j<n; j++) {        // n = m
            if (i==j) {
                // no self field
                continue;
            }
            // Calculate r'
            double rx = centx[j] - centx[i];
            double ry = centy[j] - centy[i]; 
            double rz = centz[j] - centz[i]; 
            double rmag = 1/sqrt(rx*rx + ry*ry + rz*rz);

            // Calculate cross-product 
            double jxrpx = Jy[i]*rz - Jz[i]*ry; 
            double jxrpy = Jz[i]*rx - Jx[i]*rz; 
            double jxrpz = Jx[i]*ry - Jy[i]*rx;

            // Compute contribution to field 
            Bx[j] += MU04PI * vol[i] * jxrpx * pow(rmag,3);
            By[j] += MU04PI * vol[i] * jxrpy * pow(rmag,3);
            Bz[j] += MU04PI * vol[i] * jxrpz * pow(rmag,3);
        }
    }

    return 0;
}


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
) {
    double point_creation_time, tree_build_time, moment_calc_time, bfield_calc_time;
    double total_time;
    clock_t start, end; 

    // create the octree
    // make the points first
    start = clock();
    Point **points = points_from_elements(centx, centy, centz, vol, Jx, Jy, Jz, n);
    Node *root = root_from_coords(centx, centy, centz, n);
    end = clock(); 
    point_creation_time = (double)(end-start)/CLOCKS_PER_SEC;

    // build the tree 
    start = clock();
    int success = add_points(root, points, n, 0, n);
    end = clock(); 
    tree_build_time = (double)(end-start)/CLOCKS_PER_SEC;

    // Calculate the volume*J moments 
    start = clock();
    success += calculate_moments(root);
    end = clock();
    moment_calc_time = (double)(end - start)/CLOCKS_PER_SEC;

    // recursively traverse the tree and calculate the magnetic flux density
    start = clock();
    for (size_t i=0; i<m; i++) {
        success += bfield_node_contribution(root, x[i], y[i], z[i], Bx, By, Bz, i, phi);
    }
    end = clock(); 
    bfield_calc_time = (double)(end - start)/CLOCKS_PER_SEC;

    total_time = point_creation_time + tree_build_time + moment_calc_time + bfield_calc_time;

    printf("Octree bfield calculation total time: %f s\n", total_time); 
    printf("Point creation: %.2f%%\n", 100*point_creation_time/total_time);
    printf("Tree build:     %.2f%%\n", 100*tree_build_time/total_time);
    printf("Moment calc:    %.2f%%\n", 100*moment_calc_time/total_time);
    printf("Bfield calc:    %.2f%%\n", 100*bfield_calc_time/total_time);
    return success;
}