// #include "../include/internal/csv.h"
#include "../include/thor.h"
#include <stddef.h>
#include <time.h>

const char* input_file = "data/solenoid.csv";

int main() {
    struct CsvData data = read_csv(input_file, 7, 1000, 0, ',');

    double *centx = (double*)calloc(data.nrows, sizeof(double));
    double *centy = (double*)calloc(data.nrows, sizeof(double));
    double *centz = (double*)calloc(data.nrows, sizeof(double));
    double *vol = (double*)calloc(data.nrows, sizeof(double));
    double *jx = (double*)calloc(data.nrows, sizeof(double));
    double *jy = (double*)calloc(data.nrows, sizeof(double));
    double *jz = (double*)calloc(data.nrows, sizeof(double));

    size_t n_sources = data.nrows;
    size_t ncols = data.ncols;
    size_t n_targets = 1000;
    double phi = 0.5;
    size_t i = 0; //n_targets - 1;
    Axis axis = X; 
    double line_start = 0; 
    double line_end = 1000;
    for (int i=0; i<n_sources; i++) {
        centx[i] = data.values[i*ncols+0];
        centy[i] = data.values[i*ncols+1];
        centz[i] = data.values[i*ncols+2];
        vol[i] = data.values[i*ncols+3];
        jx[i] = data.values[i*ncols+4];
        jy[i] = data.values[i*ncols+5];
        jz[i] = data.values[i*ncols+6];
    }

    printf("Solenoid problem benchmark: %li x %li\n", n_sources, n_targets);
    Line *line = new_line(axis, line_start, line_end, n_targets);
    clock_t start = clock(); 
    bfield_direct_simd(centx, centy, centz, vol, jx, jy, jz, n_sources, line->x, line->y, line->z, n_targets, line->Bx, line->By, line-> Bz, 1);
    clock_t end = clock(); 
    double elapsed_direct = (double)(end - start)/CLOCKS_PER_SEC;
    printf("Direct sum time: %.3f sec\n", elapsed_direct);
    Line *line2 = new_line(axis, line_start, line_end, n_targets);

    start = clock(); 
    bfield_octree(centx, centy, centz, vol, jx, jy, jz, n_sources, line->x, line->y, line->z, n_targets, line2->Bx, line2->By, line2-> Bz, 1, phi);
    end = clock(); 
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC; 
    printf("Speedup: %.2fx\n", elapsed_direct/elapsed);

    // Write data to file 
    
    // CsvData direct_result = 

    
    double x = line->x[i]; 
    double y = line->y[i];
    double z = line->z[i];
    double Bz_direct = line->Bz[i]; 
    double Bz_octree = line2->Bz[i];
    double err = (Bz_octree - Bz_direct) / Bz_direct;
    printf("Bz direct @ (%.3f, %.3f, %.3f) = %.6f T\n", x, y, z, Bz_direct);
    printf("Bz octree @ (%.3f, %.3f, %.3f) = %.6f T\n", x, y, z, Bz_octree);
    printf("Error: %.3f %%\n", 100*err);
}