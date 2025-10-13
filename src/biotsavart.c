#include <stdlib.h>
#include <omp.h>

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
    #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1) reduction(+:Bx[:m], By[:m], Bz[:m])
    for (size_t i=0; i<n; i++) {

        // Inner loop over obs pts 
        for (size_t j=0; j<m; j++) {
            // Calculate r'
            double rx = x[j] - centx[i];
            double ry = y[j] - centy[i]; 
            double rz = z[j] - centz[i]; 
            double rmag = 1/sqrt(rx*rx + ry*ry + rz*rz);

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
            Bx[j] += MU04PI * vol[i] * jxrpx * rmag;
            By[j] += MU04PI * vol[i] * jxrpy * rmag;
            Bz[j] += MU04PI * vol[i] * jxrpz * rmag;
        }
    }

    return 0;
}

// Naive integration of the Biot Savart law for point sources and self_fields using octree method
// B = mu0/4pi * vol * J x r' / |r'|^3
int bfield_self_naive_octree(
    const double *restrict centx, const double *restrict centy, const double *restrict centz, 
    const double *restrict vol, const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n,
    double *restrict Bx, double *restrict By, double *restrict Bz
) {

    // create the octree
    Point **points = points_from_elements(centx, centy, centz, vol, Jx, Jy, Jz, n);
    Node *root = root_from_coords(centx, centy, centz, n);
    int success = add_points(root, points, n, 0, n);

    // recursively traverse the tree and calculate the magnetic flux density
    

}