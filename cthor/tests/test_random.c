// Test a random distribution of sources 

#include "../include/thor.h"
#include <stdlib.h>
#include <time.h>

const double J = 1e8;

// Call this program like
// > test_random NUM_POINTS PHI
int main(int argc, char *argv[]) {

    // ---
    // Setup command line arguments
    // ---

    uint32_t n;
    double phi;

    if (argc == 3) {
        n = atoi(argv[1]);
        phi = (double)atof(argv[2]);
    }
    else {
        printf("Error. Incorrect number of arguments specified.\n");
        return -1;
    }

    // --- 
    // Allocate memory for sources/targets/fields
    // 
    double *x = calloc(n, sizeof(double));
    double *y = calloc(n, sizeof(double));
    double *z = calloc(n, sizeof(double));
    double *vol = calloc(n, sizeof(double));
    double *Jx = calloc(n, sizeof(double));
    double *Jy = calloc(n, sizeof(double));
    double *Jz = calloc(n, sizeof(double));
    for (int i=0; i<n; i++) {
        x[i] = drand(0.0, 1.0);
        y[i] = drand(0.0, 1.0);
        z[i] = drand(0.0, 1.0);
        vol[i] = 1e-6;
        Jx[i] = drand(-J, J);
        Jy[i] = drand(-J, J);
        Jz[i] = drand(-J, J);
    }
    double *Bx_octree = calloc(n, sizeof(double));
    double *By_octree = calloc(n, sizeof(double));
    double *Bz_octree = calloc(n, sizeof(double));
    double *Bx_direct = calloc(n, sizeof(double));
    double *By_direct = calloc(n, sizeof(double));
    double *Bz_direct = calloc(n, sizeof(double));

    // Solve 


    clock_t start = clock();
    bfield_octree(
        x, y, z, vol, Jy, Jy, Jz, n, 
        x, y, z, n, Bx_octree, By_octree, Bz_octree, 
        1, phi);
        time_t end = clock();
    double time_octree = (double)(end-start)/CLOCKS_PER_SEC;
    printf("\n");
    printf("Octree calculation time: %f s\n", time_octree);

    start = clock();
    bfield_direct_simd(
        x, y, z, vol, Jy, Jy, Jz, n, 
        x, y, z, n, Bx_direct, By_direct, Bz_direct, 
        1);

    end = clock();
    double time_direct = (double)(end-start)/CLOCKS_PER_SEC;
    printf("\n");
    printf("Direct calculation time: %f s\n", time_direct);
    printf("Speedup: %.3fx\n", time_direct / time_octree); 



    double octree_direct_error_Bx = rms_error(Bx_octree, Bx_direct, n);
    double octree_direct_error_By = rms_error(By_octree, By_direct, n);
    double octree_direct_error_Bz = rms_error(Bz_octree, Bz_direct, n);

    // printf("RMS errors relative to analytical:\n");
    printf("\n");
    printf("RMS error relative to direct:\n");
    printf("Bx octree error: %.3f %%\n", 100*octree_direct_error_Bx);
    printf("By octree error: %.3f %%\n", 100*octree_direct_error_By);
    printf("Bz octree error: %.3f %%\n", 100*octree_direct_error_Bz);

    
}