// Test a solid solenoid with the same parameters as the basic solenoid test

#include "../include/thor.h"
#include "time.h"

const double inner_radius = 0.1; 
const double outer_radius = 0.2; 
const double length = 200.0; 
const double element_size = 0.08; 
const double zcentroid = 0.0; 
const double current = 1e3; 
// const uint32_t n_targets = 1000;
const double phi = 1e-2;
const double zspan = length/100.0;

int main() {
    DenseSolenoid *ds = new_dense_solenoid(inner_radius, outer_radius, 
        length, zcentroid, current, element_size);
    uint32_t n_targets = ds->n;
    printf("n = %i\n", ds->n);

    double Bz_analytical = bfield_ideal_solenoid(ds->turns_per_unit_length, ds->current); 

    Line *line_octree = new_line(Z, -zspan/2.0, zspan/2.0, n_targets);
    Line *line_direct = new_line(Z, -zspan/2.0, zspan/2.0, n_targets);  

    time_t start = clock();
    bfield_octree(ds->x, ds->y, ds->z, ds->vol, ds->Jx, ds->Jy, ds->Jz, ds->n, 
        ds->x, ds->y, ds->z, line_octree->n, 
        line_octree->Bx, line_octree->By, line_octree->Bz, 1, phi);
    time_t end = clock();
    double time_octree = (double)(end-start)/CLOCKS_PER_SEC;
    printf("\n");
    printf("Octree calculation time: %f s\n", time_octree);

    start = clock();
    bfield_direct_simd(ds->x, ds->y, ds->z, ds->vol, ds->Jx, ds->Jy, ds->Jz, ds->n, 
        ds->x, ds->y, ds->z, line_direct->n, 
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

    double octree_direct_error_Bx = rms_error(line_octree->Bx, line_direct->Bx, line_octree->n);
    double octree_direct_error_By = rms_error(line_octree->By, line_direct->By, line_octree->n);
    double octree_direct_error_Bz = rms_error(line_octree->Bz, line_direct->Bz, line_octree->n);

    // printf("RMS errors relative to analytical:\n");
    // printf("Bz  direct error: %.3f %%\n", 100*Bz_direct_error);
    // printf("Bz  octree error: %.3f %%\n", 100*Bz_octree_error);
    printf("\n");
    printf("RMS error relative to direct:\n");
    printf("Bx octree error: %.3f %%\n", 100*octree_direct_error_Bx);
    printf("By octree error: %.3f %%\n", 100*octree_direct_error_By);
    printf("Bz octree error: %.3f %%\n", 100*octree_direct_error_Bz);
    // printf("Max field direct: (%.6f, %.6f, %.6f)\n", max(line_direct->Bx, n_targets), max(line_direct->By, n_targets), max(line_direct->Bz, n_targets));
    // printf("Min field direct: (%.6f, %.6f, %.6f)\n", min(line_direct->Bx, n_targets), min(line_direct->By, n_targets), min(line_direct->Bz, n_targets));
    // printf("Max field octree: (%.6f, %.6f, %.6f)\n", max(line_octree->Bx, n_targets), max(line_octree->By, n_targets), max(line_octree->Bz, n_targets));
    // printf("Min field octree: (%.6f, %.6f, %.6f)\n", min(line_octree->Bx, n_targets), min(line_octree->By, n_targets), min(line_octree->Bz, n_targets));
    
    free_dense_solenoid(ds);
    free_line(line_octree);
    free_line(line_direct);
}