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
double zspan = 5.0;            // Length of line to calculate results on
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
    // Create target points along the z-axis and allocate results arrays
    // --- 

    // Compress points to region where analytical calculation is accurate
    double zcenter = zspan / 2.0; 
    double line_start = zcenter - ( zspan / (2.0 * compression) );
    double line_end = line_start + zspan / compression;

    Line *line_direct = new_line(Z, line_start, line_end, n_targets);
    Line *line_octree = new_line(Z, line_start, line_end, n_targets);
    double *Bz_analytical = calloc(n_targets, sizeof(double));


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
    double zstep = 0.0;

    int elements_per_turn = 100;
    side_length = zspan / ((double)n_sources / (double)elements_per_turn);
    area = side_length*side_length;
    v = side_length*area;        
    theta_step = 2.0*PI/(double)(elements_per_turn - 1);
    zstep = zspan/(double)n_sources; 
    J = I / area;
    

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
    // Should not change along the axis
    double Bz_solenoid = bfield_ideal_solenoid(100.0, I); 
    for (size_t i=0; i<n_targets; i++) {
        Bz_analytical[i] = Bz_solenoid;
    }

    // Octree calculation (timed internally but repeated here for calculation of speedup)
    time_t start = clock();
    bfield_octree(
        xs, ys, zs, vol, Jx, Jy, Jz, n_sources, 
        line_octree->x, line_octree->y, line_octree->z, n_targets, 
        line_octree->Bx, line_octree->By, line_octree->Bz, 1, phi);
    time_t end = clock();
    double time_octree = (double)(end-start)/CLOCKS_PER_SEC;
    printf("\n");
    printf("Octree calculation time: %f s\n", time_octree);

    // Direct summation of all points
    start = clock();
    bfield_direct(
        xs, ys, zs, vol, Jx, Jy, Jz, n_sources, 
        line_direct->x, line_direct->y, line_direct->z, n_targets, 
        line_direct->Bx, line_direct->By, line_direct->Bz, 1);
    end = clock(); 
    double time_direct = (double)(end-start)/CLOCKS_PER_SEC;
    printf("Direct calculation time: %f s\n", time_direct);
    printf("Speedup: %.3fx\n", time_direct/time_octree);
    printf("\n");
    
    // --- 
    // Compare solutions 
    // --- 

    double Bz_direct_error = rms_error(line_direct->Bz, Bz_analytical, n_targets); 
    double Bz_octree_error = rms_error(line_octree->Bz, Bz_analytical, n_targets);
    double octree_direct_error = rms_error(line_octree->Bz, line_direct->Bz, n_targets);

    printf("RMS errors relative to analytical:\n");
    printf("Bz  direct error: %.3f %%\n", 100*Bz_direct_error);
    printf("Bz  octree error: %.3f %%\n", 100*Bz_octree_error);
    printf("\n");
    printf("RMS error relative to direct:\n");
    printf("Bz octree error: %.3f %%\n", 100*octree_direct_error);

    int j = n_targets/2;
    printf("Analytical field at z = %.3f: %.3f T\n", line_direct->z[j], Bz_analytical[j]);
    printf("Direct     field at z = %.3f: %.3f T\n", line_direct->z[j], line_direct->z[j]);
    printf("Octree     field at z = %.3f: %.3f T\n", line_direct->z[j], line_octree->z[j]);


    // Ctrl+F gives 8 `calloc()` and 8 `free()` commands in this function
    free(Bz_analytical); 
    free_line(line_direct);
    free_line(line_octree);
    free(xs); free(ys); free(zs); free(vol); free(Jx); free(Jy); free(Jz);
}