// Confirm that the code accurately calculates the magnetic field at the 
// center of a current carrying wire loop. The loop is in the xy plane. 
// The observation points are on a line on the z-axis. The magnetic field
// should be only in the z-direction.

#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include "../src/biotsavart.h"
#include "../src/utils.h"


// Constant parameters for the problem
// (these are not that important to change between tests)
const double I = 1e6;           // Total current in loop 
const double R = 1.0;           // Radius of loop
const double z_span = 10.0;     // Length of line to calculate results on


// call this function like:
//  ./test_loop.o NUM_SOURCES NUM_TARGETS NUM_THREADS PHI
int main(int argc, char *argv[]) {

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

    // Create target points along the z-axis, centered about the origin
    const double *xt = calloc(n_targets, sizeof(double));
    const double *yt = calloc(n_targets, sizeof(double));
    double *zt = calloc(n_targets, sizeof(double));
    double zstep = z_span / (double)(n_targets - 1);
    zt[0] = -z_span/2.0;
    for (uint32_t i=1; i<n_targets; i++) {
        zt[i] = zt[i-1] + zstep;
    }

    // Allocate results arrays
    double *Bz_analytical = calloc(n_targets, sizeof(double));
    double *Bx_direct = calloc(n_targets, sizeof(double));
    double *By_direct = calloc(n_targets, sizeof(double));
    double *Bz_direct = calloc(n_targets, sizeof(double));
    double *Bx_octree = calloc(n_targets, sizeof(double));
    double *By_octree = calloc(n_targets, sizeof(double));
    double *Bz_octree = calloc(n_targets, sizeof(double));

    // Create source points 
    double theta_step = 2*PI / (double)(n_sources - 1); 
    double step = theta_step*R;
    double side_length = R / 1000; 
    double area = side_length * side_length;
    double total_volume = area * 2*PI*R;
    double v = total_volume/n_sources;
    double J = I / area;
    double *xs = calloc(n_sources, sizeof(double));
    double *ys = calloc(n_sources, sizeof(double));
    double *zs = calloc(n_sources, sizeof(double));
    double *vol = calloc(n_sources, sizeof(double));
    double *Jx = calloc(n_sources, sizeof(double));
    double *Jy = calloc(n_sources, sizeof(double));
    double *Jz = calloc(n_sources, sizeof(double));
    double theta = 0; 
    for (size_t i = 0; i<n_sources; i++) {
        xs[i] = R*cos(theta); 
        ys[i] = R*sin(theta);
        vol[i] = v;
        Jx[i] = -J*sin(theta); 
        Jy[i] = J*cos(theta);
        theta += theta_step;
    }

    // Compute solutions

    // Analytical solution to test both direct and octree methods against
    bfield_loop_axis(zt, n_targets, I, R, Bz_analytical); 

    // Direct summation of all points
    time_t start = clock();
    bfield_direct(
        xs, ys, zs, vol, Jx, Jy, Jz, n_sources, xt, yt, zt, n_targets, 
        Bx_direct, By_direct, Bz_direct, 1);
    time_t end = clock(); 

    printf("Direct calculation time: %f s\n", (double)(end-start)/CLOCKS_PER_SEC);

    // Octree calculation (timed internally)
    bfield_octree(xs, ys, zs, vol, Jx, Jy, Jz, n_sources, xt, yt, zt, n_targets, 
        Bx_octree, By_octree, Bz_octree, 1, phi);
    
    // Compare solutions 
    double Bxy_direct_error = 0.0;      
    double Bxy_octree_error = 0.0; 
    double Bz_direct_error = rms_error(Bz_direct, Bz_analytical, n_targets); 
    double Bz_octree_error = rms_error(Bz_octree, Bz_analytical, n_targets);

    printf("Average errors:\n");
    printf("Bxy direct error: %.3f %%\n", 100*Bxy_direct_error);
    printf("Bxy octree error: %.3f %%\n", 100*Bxy_octree_error);
    printf("Bz  direct error: %.3f %%\n", 100*Bz_direct_error);
    printf("Bz  octree error: %.3f %%\n", 100*Bz_octree_error);

    int j = 0;
    printf("Analytical field at z = %.3f: %.3f T\n", zt[j], Bz_analytical[j]);
    printf("Direct     field at z = %.3f: %.3f T\n", zt[j], Bz_direct[j]);
    printf("Octree     field at z = %.3f: %.3f T\n", zt[j], Bz_octree[j]);


    free(xt); free(yt); free(zt); 
    free(Bz_analytical); 
    free(Bx_direct); free(By_direct); free(Bz_direct); 
    free(Bx_octree); free(By_octree); free(Bz_octree);
    free(xs); free(ys); free(zs); free(vol); free(Jx); free(Jy); free(Jz);
}