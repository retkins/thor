// Test a solid solenoid with the same parameters as the basic solenoid test

#include "../include/thor.h"
#include "time.h"

// Test like its a circular ring
const double inner_radius = 1.0; 
const double outer_radius = 1.1; 
const double length = 1.0; 
const double element_size = 0.015; 
const double zcentroid = 0.0; 
const double current = 1e3; 

const double phi = 1e-1;
const double zspan = length/100.0;

int main() {
    DenseSolenoid *ds = new_dense_solenoid(inner_radius, outer_radius, 
        length, zcentroid, current, element_size);
    const uint32_t n_targets = ds->n;
    printf("n = %i\n", ds->n);
    printf("m = %i\n", n_targets);

    double Bz_analytical = bfield_ideal_solenoid(ds->turns_per_unit_length, ds->current); 
    double R = (outer_radius - inner_radius)/2.0;
    // double Bz_analytical = bfield_loop_axis(line_octree->z, n_targets, current, R, double *restrict Bz)
    Line *line_octree = new_line(Z, -zspan/2.0, zspan/2.0, n_targets);
    Line *line_direct = new_line(Z, -zspan/2.0, zspan/2.0, n_targets);  

    time_t start = clock();
    bfield_octree(ds->x, ds->y, ds->z, ds->vol, ds->Jx, ds->Jy, ds->Jz, ds->n, 
        line_octree->x, line_octree->y, line_octree->z, line_octree->n, 
        line_octree->Bx, line_octree->By, line_octree->Bz, 1, phi);
    time_t end = clock();
    double time_octree = (double)(end-start)/CLOCKS_PER_SEC;
    printf("\n");
    printf("Octree calculation time: %f s\n", time_octree);

    start = clock();
    bfield_direct_simd(ds->x, ds->y, ds->z, ds->vol, ds->Jx, ds->Jy, ds->Jz, ds->n, 
        line_direct->x, line_direct->y, line_direct->z, line_direct->n, 
        line_direct->Bx, line_direct->By, line_direct->Bz, 1);
    end = clock();
    double time_direct = (double)(end-start)/CLOCKS_PER_SEC;
    printf("\n");
    printf("Direct calculation time: %f s\n", time_direct);
    printf("Speedup: %.3fx\n", time_direct / time_octree); 

    int j = n_targets/2;
    printf("Bz analytical: %.6f T\n", Bz_analytical);
    printf("Bz octree:     %.6f T\n", line_octree->Bz[j]);
    printf("Bz direct:     %.6f T\n", line_direct->Bz[j]);

    double octree_direct_error = rms_error(line_octree->Bz, line_direct->Bz, line_octree->n);

    // printf("RMS errors relative to analytical:\n");
    // printf("Bz  direct error: %.3f %%\n", 100*Bz_direct_error);
    // printf("Bz  octree error: %.3f %%\n", 100*Bz_octree_error);
    printf("\n");
    printf("RMS error relative to direct:\n");
    printf("Bz octree error: %.3f %%\n", 100*octree_direct_error);

    free_dense_solenoid(ds);
    free_line(line_octree);
    free_line(line_direct);
}