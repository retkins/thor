#include <stdlib.h>
#include <time.h>
#include <math.h>
// #include <omp.h>

#include "../include/internal/biotsavart.h"
#include "../include/internal/octree.h"
#include "../include/internal/utils.h"
#include "../include/internal/simd.h"

// --- 
// Analytical solutions 
// --- 


double bfield_ideal_solenoid(double turns_per_unit_length, double current) {
    return MU0 * turns_per_unit_length * current; 
}


int bfield_loop_axis(
    const double *restrict z, size_t n, 
    const double I, const double R, 
    double *restrict Bz
) {
    double R2 = R*R;
    for (size_t i=0; i<n; i++) {
        Bz[i] = MU0 * I * R2 / (2 * pow(z[i]*z[i] + R2, 1.5));
    }

    return 0;
}


// --- 
// Discrete integration solutions 
// ---


int bfield_direct(
    const double *restrict cx, const double *restrict cy, const double *restrict cz, 
    const double *restrict vol, const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n_sources,
    const double *restrict x, const double *restrict y, const double *restrict z, 
    size_t n_targets,
    double *restrict Bx, double *restrict By, double *restrict Bz, 
    int nthreads
) {

    const double ONE_OVER_FOUR_THIRDS_PI = 1 / (4.0 * PI / 3.0 );
    const double ONE_THIRD = 1.0/3.0;

    // Outer loop over sources
    // #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1) reduction(+:Bx[:m], By[:m], Bz[:m])
    for (size_t i=0; i<n_sources; i++) {

        // volume = 4/3 * pi * R^3
        double R = powf(ONE_OVER_FOUR_THIRDS_PI * vol[i], ONE_THIRD);

        // Inner loop over target pts 
        for (size_t j=0; j<n_targets; j++) {
            // Calculate r'
            double rx = x[j] - cx[i];
            double ry = y[j] - cy[i]; 
            double rz = z[j] - cz[i]; 
            double rmag = sqrt(rx*rx + ry*ry + rz*rz);
            double rmag3 = rmag*rmag*rmag;
            double inv_rmag3;
            if (rmag > R) { inv_rmag3 = 1 / rmag3; }
            else { inv_rmag3 = powf(rmag / R, 3.0); }

            // Calculate cross-product 
            double jxrpx = Jy[i]*rz - Jz[i]*ry; 
            double jxrpy = Jz[i]*rx - Jx[i]*rz; 
            double jxrpz = Jx[i]*ry - Jy[i]*rx;

            // Compute contribution to field 
            Bx[j] += MU04PI * vol[i] * jxrpx * inv_rmag3;
            By[j] += MU04PI * vol[i] * jxrpy * inv_rmag3;
            Bz[j] += MU04PI * vol[i] * jxrpz * inv_rmag3;
        }
    }

    return 0;
}


int bfield_direct_simd(
    const double *restrict cx, const double *restrict cy, const double *restrict cz, 
    const double *restrict vol, const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n_sources,
    const double *restrict x, const double *restrict y, const double *restrict z, 
    size_t n_targets,
    double *restrict Bx, double *restrict By, double *restrict Bz, 
    int nthreads
) {

    const double ONE_OVER_FOUR_THIRDS_PI = 1 / (4.0 * PI / 3.0 );
    const double ONE_THIRD = 1.0/3.0;
    vpd _MU04PI = setpd(MU04PI);

    // Outer loop over sources
    // #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1) reduction(+:Bx[:m], By[:m], Bz[:m])
    for (size_t i=0; i<n_sources; i++) {

        // volume = 4/3 * pi * R^3
        double R = powf(ONE_OVER_FOUR_THIRDS_PI * vol[i], ONE_THIRD);
        vpd _R = setpd(R);
        vpd R3 = mulpd(mulpd(_R, _R), _R);

        // Vector loads before loop 
        vpd _Jx = setpd(Jx[i]); 
        vpd _Jy = setpd(Jy[i]); 
        vpd _Jz = setpd(Jz[i]);
        vpd _MU04PI_vol = mulpd(_MU04PI, setpd(vol[i]));

        // Inner loop over target pts 
        // SIMD strides first
        size_t j=0; 
        for (; j<n_targets; j+=VLEND) {

            // calculate r' = t[j] - s[i]
            vpd rx = subpd(loadpd(&x[j]), setpd(cx[i]));
            vpd ry = subpd(loadpd(&y[j]), setpd(cy[i]));
            vpd rz = subpd(loadpd(&z[j]), setpd(cz[i]));
            vpd rmag = mag3pd(rx, ry, rz);
            vpd rmag3 = mulpd(mulpd(rmag, rmag), rmag);

            // Branchless means of reducing the field inside the element
            // (graceful singularity handling)
            vpd outside = invpd(rmag3); 
            vpd inside = mulpd(rmag3, invpd(R3));
            vpdu mask = cmpgtpd(rmag, _R); 
            vpd inv_rmag3 = blendpd(outside, inside, mask);

            // Cross-product 
            vpd jxrpx = subpd(mulpd(_Jy,rz), mulpd(_Jz,ry));
            vpd jxrpy = subpd(mulpd(_Jz,rx), mulpd(_Jx,rz));
            vpd jxrpz = subpd(mulpd(_Jx,ry), mulpd(_Jy,rx));

            // Compute bfield contribution
            vpd _Bx = mulpd(mulpd(_MU04PI_vol, jxrpx), inv_rmag3);
            vpd _By = mulpd(mulpd(_MU04PI_vol, jxrpy), inv_rmag3);
            vpd _Bz = mulpd(mulpd(_MU04PI_vol, jxrpz), inv_rmag3);

            vpd _Bx_old = loadpd(&Bx[j]);
            vpd _By_old = loadpd(&By[j]);
            vpd _Bz_old = loadpd(&Bz[j]);

            // Store 
            storepd(&Bx[j], addpd(_Bx, _Bx_old));
            storepd(&By[j], addpd(_By, _By_old));
            storepd(&Bz[j], addpd(_Bz, _Bz_old));
        }
        
        // Scalar fallback to handle remainder
        for (; j<n_targets; j++) {
            // Calculate r'
            double rx = x[j] - cx[i];
            double ry = y[j] - cy[i]; 
            double rz = z[j] - cz[i]; 
            double rmag = sqrt(rx*rx + ry*ry + rz*rz);
            double rmag3 = rmag*rmag*rmag;
            double inv_rmag3;
            if (rmag > R) { inv_rmag3 = 1 / rmag3; }
            else { inv_rmag3 = powf(rmag / R, 3.0); }

            // Calculate cross-product 
            double jxrpx = Jy[i]*rz - Jz[i]*ry; 
            double jxrpy = Jz[i]*rx - Jx[i]*rz; 
            double jxrpz = Jx[i]*ry - Jy[i]*rx;

            // Compute contribution to field 
            Bx[j] += MU04PI * vol[i] * jxrpx * inv_rmag3;
            By[j] += MU04PI * vol[i] * jxrpy * inv_rmag3;
            Bz[j] += MU04PI * vol[i] * jxrpz * inv_rmag3;
        }
    }

    return 0;
}


int bfield_direct_self(
    const double *restrict cx, const double *restrict cy, const double *restrict cz, 
    const double *restrict vol, const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n, // sources and target points
    double *restrict Bx, double *restrict By, double *restrict Bz,
    int nthreads
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
            double rx = cx[j] - cx[i];
            double ry = cy[j] - cy[i]; 
            double rz = cz[j] - cz[i]; 
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


int bfield_octree(
    const double *restrict cx, const double *restrict cy, const double *restrict cz, 
    const double *restrict vol, const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n_sources, 
    const double *restrict x, const double *restrict y, const double *restrict z, 
    size_t n_targets,
    double *restrict Bx, double *restrict By, double *restrict Bz, 
    int nthreads, 
    double phi
) {
    // Timing variables (helpful for performance engineering)
    double point_creation_time, tree_build_time, moment_calc_time, bfield_calc_time;
    double total_time;
    clock_t start, end; 

    // --- 
    // Create the octree 
    // ---

    start = clock();

    // TODO: there's a magic number for allocation of nodes here; need to 
    // improve the reallocation behavior of the nodes array

    Point *points = points_from_elements(cx, cy, cz, vol, Jx, Jy, Jz, n_sources);
    NodeAllocation *nodes = allocate_nodes(10*n_sources);        
    Node *root = root_from_coords(nodes, cx, cy, cz, n_sources);
    
    end = clock(); 
    point_creation_time = (double)(end-start)/CLOCKS_PER_SEC;
    // It's interesting that the time to create the points is usually slightly
    // more than the upwards pass (moment calc)

    // Build the tree by recursively placing points within it
    start = clock();
    int success = add_points(nodes, root, points, n_sources, 0, n_sources);
    end = clock(); 
    tree_build_time = (double)(end-start)/CLOCKS_PER_SEC;

    // Calculate the volume*J moments (upwards pass)
    start = clock();
    success += calculate_moments(root);
    end = clock();
    moment_calc_time = (double)(end - start)/CLOCKS_PER_SEC;

    // recursively traverse the tree and calculate the magnetic flux density
    // at each target point (downwards pass)
    start = clock();
    for (size_t i=0; i<n_targets; i++) {
        success += bfield_node_contribution(root, x[i], y[i], z[i], Bx, By, Bz, i, phi);
    }
    end = clock(); 
    bfield_calc_time = (double)(end - start)/CLOCKS_PER_SEC;

    total_time = point_creation_time + tree_build_time + moment_calc_time + bfield_calc_time;

    printf("Octree bfield calculation completed using phi = %.3f\n", phi);
    printf("Octree bfield calculation total time: %f s\n", total_time); 
    printf("Point creation: %.3f s (%.2f%%)\n", point_creation_time, 100*point_creation_time/total_time);
    printf("Tree build:     %.3f s (%.2f%%)\n", tree_build_time, 100*tree_build_time/total_time);
    printf("Moment calc:    %.3f s (%.2f%%)\n", moment_calc_time, 100*moment_calc_time/total_time);
    printf("Bfield calc:    %.3f s (%.2f%%)\n", bfield_calc_time, 100*bfield_calc_time/total_time);
    
    deallocate_nodes(nodes); free(points);
    return success;
}


// Naive O(N^2) integration of the Biot Savart Law for finite-length filaments
// 
// Reference:
// https://freestatelabs.github.io/Wired.jl/dev/theory/
int bfield_wire(
    const double *restrict cx, const double *restrict cy, const double *restrict cz, 
    const double *restrict vol, 
    const double *restrict Jx, const double *restrict Jy, const double *restrict Jz, 
    size_t n_sources, 
    const double *restrict x, const double *restrict y, const double *restrict z, 
    double *restrict Bx, double *restrict By, double *restrict Bz, 
    size_t n_targets
) {
    for (size_t i=0; i<n_sources; i++) {
        // Compute the start and end points 
        double Jmag = sqrt(Jx[i]*Jx[i] + Jy[i]*Jy[i] + Jz[i]*Jz[i]); 
        double uJx=0; double uJy=0; double uJz=0;
        if (Jmag > 1e-8 ) {
            uJx = Jx[i]/Jmag; uJy = Jy[i]/Jmag; uJz = Jz[i]/Jmag;
        }

        // Assume element shape is a cube: vol = 8*R^3, or vol = length^3
        double length = cbrt(vol[i]);
        double R = length/2.0; 
        double ax0 = cx[i] - R*uJx; 
        double ax1 = ax0 + length*uJx;
        double ay0 = cy[i] - R*uJy; 
        double ay1 = ay0 + length*uJy;
        double az0 = cz[i] - R*uJz; 
        double az1 = az0 + length*uJz;
        double ax = ax1 - ax0; double ay = ay1 - ay0; double az = az1 - az0; 
        // Compute total current 
        double area = length*length; 
        double I_mu0_4pi = (area*Jmag) * MU04PI;
        
        // Loop over all targets
        for (size_t j=0; j<n_targets; j++) {

            // Vectors 
            double bx = ax0 - x[j]; 
            double by = ay0 - y[j]; 
            double bz = az0 - z[j]; 
            double cx = ax1 - x[j]; 
            double cy = ay1 - y[j]; 
            double cz = az1 - z[j]; 

            // Dot products and magnitudes 
            double a_dot_c = ax*cx + ay*cy + az*cz; 
            double a_dot_b = ax*bx + ay*by + az*bz;
            double c_mag = sqrt(cx*cx + cy*cy + cz*cz); 
            double b_mag = sqrt(bx*bx + by*by + bz*bz);

            // Cross products 
            double c_cross_a_x = cy*az - cz*ay; 
            double c_cross_a_y = cz*ax - cx*az; 
            double c_cross_a_z = cx*ay - cy*ax; 
            double c_cross_a_mag2 = c_cross_a_x*c_cross_a_x + c_cross_a_y*c_cross_a_y + c_cross_a_z*c_cross_a_z;

            // Contributions to field
            double coeff = I_mu0_4pi * ((a_dot_c/c_mag) - (a_dot_b/b_mag)) / c_cross_a_mag2;
            Bx[j] += coeff*c_cross_a_x;
            By[j] += coeff*c_cross_a_y; 
            Bz[j] += coeff*c_cross_a_z;
        }
    }

    return 0;
}