/*  Confirm that the code accurately calculates the magnetic field at the 
    center of a current carrying wire loop. The loop is in the xy plane. 
    The observation points are on a line on the z-axis. The magnetic field
    should be only in the z-direction.

    Note: this test can be performed for either a loop or a solenoid. 
    TODO: make the solenoid a separate test. 
*/

#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include "../include/thor.h"

// Constant parameters for the problem
const double I = 1e3;           // Total current in loop 
const double R = 1.0/(2*PI);    // Radius of loop
double z_span = 5.0;            // Length of line to calculate results on
bool solenoid = true;
double compression = 2.0;       // 'compress' target pts at center of solenoid
// if the solenoid is requested, then the element size changes to 1 cm 
// and 100 elements/turn


// call this function like:
//  ./test_loop.o NUM_SOURCES NUM_TARGETS NUM_THREADS PHI
int main(int argc, char *argv[]) {

    // ---
    // Setup command line arguments
    // ---

    uint32_t n_sources, n_targets, n_threads;
    double phi;

    if (argc == 5) {
        n_sources = atoi(argv[1]);
        n_targets = atoi(argv[2]);
        n_threads = atoi(argv[3]);
        phi = (double)atof(argv[4]);
    }
    else {
        printf("Error. Incorrect number of arguments specified.\n");
        return -1;
    }

    // ---
    // Create target points along the z-axis
    // --- 

    double *xt = calloc(n_targets, sizeof(double));
    double *yt = calloc(n_targets, sizeof(double));
    double *zt = calloc(n_targets, sizeof(double));
    double zstep = z_span / (double)(n_targets - 1);
    zt[0] = -z_span/2.0;        // center about the origin for current loop
    if (solenoid) {
        // Compress points to region where analytical calculation is accurate
        zt[0] = 0.0;
        double zc = z_span / 2.0; 
        zt[0] = zc - z_span / (2.0 * compression);
        zstep = (z_span / compression) / (double)(n_targets - 1);
    }
    for (uint32_t i=1; i<n_targets; i++) {
        zt[i] = zt[i-1] + zstep;
    }

    // ---
    // Allocate results arrays
    // ---
    double *Bz_analytical = calloc(n_targets, sizeof(double));
    double *Bx_direct = calloc(n_targets, sizeof(double));
    double *By_direct = calloc(n_targets, sizeof(double));
    double *Bz_direct = calloc(n_targets, sizeof(double));
    double *Bx_octree = calloc(n_targets, sizeof(double));
    double *By_octree = calloc(n_targets, sizeof(double));
    double *Bz_octree = calloc(n_targets, sizeof(double));

    // --- 
    // Create source points 
    // --- 

    double *xs = calloc(n_sources, sizeof(double));
    double *ys = calloc(n_sources, sizeof(double));
    double *zs = calloc(n_sources, sizeof(double));
    double *vol = calloc(n_sources, sizeof(double));
    double *Jx = calloc(n_sources, sizeof(double));
    double *Jy = calloc(n_sources, sizeof(double));
    double *Jz = calloc(n_sources, sizeof(double));
    // Need to do some extra work to place points in a spiral 
    double theta_step = 2*PI / (double)(n_sources - 1); 
    double step = theta_step*R;
    double side_length = R / 1000; 
    double area = side_length * side_length;
    double total_volume = area * 2*PI*R;
    double v = total_volume/n_sources;
    double J = I / area;
    double z = 0.0;
    double theta = 0; 
    zstep = 0.0;

    // Create a coil instead of a loop for testing octree times
    if (solenoid)
    {
        int elements_per_turn = 100;
        side_length = z_span / ((double)n_sources / (double)elements_per_turn);
        area = side_length*side_length;
        v = side_length*area;        
        theta_step = 2.0*PI/(double)(elements_per_turn - 1);
        zstep = z_span/(double)n_sources; 
        J = I / area;
    }

    for (size_t i = 0; i<n_sources; i++) {
        xs[i] = R*cos(theta); 
        ys[i] = R*sin(theta);
        zs[i] = z;
        vol[i] = v;
        Jx[i] = -J*sin(theta); 
        Jy[i] = J*cos(theta);
        theta += theta_step;
        z += zstep;
    }

    // ---
    // Compute solutions
    // ---

    // Analytical solution to test both direct and octree methods against
    if (solenoid) {
        // Should not change along the axis
        double Bz_solenoid = MU0 * 100.0 * I; 
        for (size_t i=0; i<n_targets; i++) {
            Bz_analytical[i] = Bz_solenoid;
        }
    }
    else { 
        bfield_loop_axis(zt, n_targets, I, R, Bz_analytical); 
    }

    // Octree calculation (timed internally but repeated here for calculation of speedup)
    time_t start = clock();
    bfield_octree(xs, ys, zs, vol, Jx, Jy, Jz, n_sources, xt, yt, zt, n_targets, 
        Bx_octree, By_octree, Bz_octree, 1, phi);
    time_t end = clock();
    double time_octree = (double)(end-start)/CLOCKS_PER_SEC;
    printf("\n");
    printf("Octree calculation time: %f s\n", time_octree);

    // Direct summation of all points
    start = clock();
    bfield_direct(
        xs, ys, zs, vol, Jx, Jy, Jz, n_sources, xt, yt, zt, n_targets, 
        Bx_direct, By_direct, Bz_direct, 1);
    end = clock(); 
    double time_direct = (double)(end-start)/CLOCKS_PER_SEC;
    printf("Direct calculation time: %f s\n", time_direct);
    printf("Speedup: %.3fx\n", time_direct/time_octree);
    printf("\n");
    
    // --- 
    // Compare solutions 
    // --- 

    double Bz_direct_error = rms_error(Bz_direct, Bz_analytical, n_targets); 
    double Bz_octree_error = rms_error(Bz_octree, Bz_analytical, n_targets);
    double octree_direct_error = rms_error(Bz_octree, Bz_direct, n_targets);

    printf("RMS errors relative to analytical:\n");
    printf("Bz  direct error: %.3f %%\n", 100*Bz_direct_error);
    printf("Bz  octree error: %.3f %%\n", 100*Bz_octree_error);
    printf("\n");
    printf("RMS error relative to direct:\n");
    printf("Bz octree error: %.3f %%\n", 100*octree_direct_error);

    int j = n_targets/2;
    printf("Analytical field at z = %.3f: %.3f T\n", zt[j], Bz_analytical[j]);
    printf("Direct     field at z = %.3f: %.3f T\n", zt[j], Bz_direct[j]);
    printf("Octree     field at z = %.3f: %.3f T\n", zt[j], Bz_octree[j]);


    // Ctrl+F gives 17 `calloc()` and 17 `free()` commands in this function
    free(xt); free(yt); free(zt); 
    free(Bz_analytical); 
    free(Bx_direct); free(By_direct); free(Bz_direct); 
    free(Bx_octree); free(By_octree); free(Bz_octree);
    free(xs); free(ys); free(zs); free(vol); free(Jx); free(Jy); free(Jz);
}