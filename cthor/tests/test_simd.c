/*  Using NEON intrinsics is about 4.7x faster than scalar clang code
*/

#include "../include/thor.h"
#include <time.h>
#include <stdio.h>

int main() {
    double radius = 1.0;
    double length = 10.0;
    double zcentroid = 0.0;
    double current = 1e3; 
    double turns_per_unit_length = 100.0;
    double n = 10000;
    double zspan = 1.0;

    Solenoid *solenoid = new_solenoid(radius, length, zcentroid, current, turns_per_unit_length, n);

    Line *line_direct = new_line(Z, -zspan/2.0, zspan/2.0, n);
    Line *line_simd = new_line(Z, -zspan/2.0, zspan/2.0, n);

    clock_t start = clock();
    bfield_direct(solenoid->x, solenoid->y, solenoid->z, solenoid->vol, 
        solenoid->Jx, solenoid->Jy, solenoid->Jz, n, 
        line_direct->x, line_direct->y, line_direct->z, n, 
        line_direct->Bx, line_direct->By, line_direct->Bz, 1);
    clock_t end = clock();

    double direct_time = (double)(end - start) / CLOCKS_PER_SEC;

    start = clock();
    bfield_direct_simd(solenoid->x, solenoid->y, solenoid->z, solenoid->vol, 
        solenoid->Jx, solenoid->Jy, solenoid->Jz, n, 
        line_simd->x, line_simd->y, line_simd->z, n, 
        line_simd->Bx, line_simd->By, line_simd->Bz, 1);
    end = clock();
    double simd_time = (double)(end - start) / CLOCKS_PER_SEC;

    double Bz_error = rms_error(line_simd->Bz, line_direct->Bz, n); 

    printf("Direct time: %.6f s\n", direct_time); 
    printf("SIMD   time: %.6f s\n", simd_time); 
    printf("Speedup: %.3fx\n", direct_time/simd_time);
    printf("Bz error:    %.6f %%\n", 100*Bz_error);

    free_solenoid(solenoid); 
    free_line(line_direct);
    free_line(line_simd);
}