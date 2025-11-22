#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h> 

#include "../include/internal/utils.h"
#include "../include/internal/problems.h"
#include "../include/internal/biotsavart.h"


// ---
// Simple current-carrying circular loop 
// ---

Loop *new_loop(double radius, double current, uint32_t n) {
    Loop *loop = calloc(1, sizeof(Loop));

    if (loop != NULL) {
        loop->radius = radius; 
        loop->current = current; 
        loop->n = n; 

        loop->x = zeros(n);
        loop->y = zeros(n);
        loop->z = zeros(n);
        loop->vol = zeros(n);
        loop->Jx = zeros(n);
        loop->Jy = zeros(n);
        loop->Jz = zeros(n);

        if (
            (loop->x != NULL) && (loop->y != NULL) && (loop->z != NULL) &&
            (loop->vol != NULL) && 
            (loop->Jx != NULL) && (loop->Jy != NULL) && (loop->Jz != NULL) 
        ) { 
            // TODO: for larger values of `n`, the element size dimishes to zero
            double theta = 0.0; 
            double theta_step = 2*PI/(double)n;
            double perimeter = 2*PI*radius; 
            double side_length = perimeter/(double)n;
            double element_area = side_length * side_length;
            double element_volume = element_area * side_length;
            double element_jdensity = current/element_area; 
            
            for (uint32_t i=1; i<n; i++) {
                loop->x[i] = radius*cos(theta); 
                loop->y[i] = radius*sin(theta);
                loop->vol[i] = element_volume;
                loop->Jx[i] = -element_jdensity*sin(theta);
                loop->Jy[i] = element_jdensity*cos(theta);
                theta += theta_step;
            }
            return loop;
        }
        else {
            // calloc() initialized all pointers to NULL, so we can free the 
            // unallocated pointers without trouble
            free_loop(loop);
            return NULL; 
        }
    
    }
    else {
            fprintf(stderr, "Error in `new_loop()`: could not allocate memory.\n");
            free_loop(loop);
            return NULL;
        }
}


void free_loop(Loop *loop) {
    if (loop != NULL) {
        free(loop->x);
        free(loop->y);
        free(loop->z);
        free(loop->vol);
        free(loop->Jx);
        free(loop->Jy);
        free(loop->Jz);
        free(loop);
        loop = NULL;
    }   
}


// ---
// Simple solenoid (one element thick in radial direction)
// ---


Solenoid *new_solenoid(
    double radius, double length, double zcentroid, double current, 
    double turns_per_unit_length, uint32_t n
) {
    Solenoid *solenoid = calloc(1, sizeof(Solenoid));

    if (solenoid != NULL) {
        solenoid->radius = radius; 
        solenoid->length = length;
        solenoid->zcentroid = zcentroid;
        solenoid->current = current;
        solenoid->turns_per_unit_length = turns_per_unit_length;
        solenoid->n = n;

        solenoid->x = zeros(n); solenoid->y = zeros(n); solenoid->z = zeros(n);
        solenoid->vol = zeros(n);
        solenoid->Jx = zeros(n); solenoid->Jy = zeros(n); solenoid->Jz = zeros(n);

        if ( 
            (solenoid->x != NULL) && (solenoid->y != NULL) && \
            (solenoid->z != NULL) && (solenoid->vol != NULL) && \
            (solenoid->Jx != NULL) && (solenoid->Jy != NULL) && \
            (solenoid->Jz != NULL) 
        ) {

            // Determine how to spiral-wind the solenoid
            double turns = turns_per_unit_length * length;      // note non-integer
            double total_theta = 2 * PI * turns;
            double theta_step = total_theta / (double)(n - 1);
            double zstep = length / (double)(n - 1);

            // Determine characteristics of the elements 
            double circumference = 2 * PI * radius;
            double total_winding_length = circumference * turns_per_unit_length * length;   // roughly
            double element_size = total_winding_length / (double)n;
            double element_area = element_size * element_size;
            double element_volume = element_area * element_size;
            double element_jdensity = current / element_area;

            // Start with first element 
            solenoid->x[0] = radius; 
            solenoid->z[0] = zcentroid - length/2.0;
            solenoid->vol[0] = element_volume;
            solenoid->Jy[0] = element_jdensity;

            double theta = theta_step;
            for (size_t i=1; i<n; i++) {
                solenoid->x[i] = radius * cos(theta); 
                solenoid->y[i] = radius * sin(theta);
                solenoid->z[i] = solenoid->z[i-1] + zstep; 
                solenoid->vol[i] = element_volume;
                solenoid->Jx[i] = -element_jdensity * sin(theta);
                solenoid->Jy[i] = element_jdensity * cos(theta); 
                theta += theta_step;
                // printf("x, y, z = %.6f, %.6f, %.6f\n", solenoid->x[i], solenoid->y[i], solenoid->z[i]);
                // printf("Jx, Jy, Jz = %.6f, %.6f, %.6f\n", solenoid->Jx[i], solenoid->Jy[i], solenoid->Jz[i]);
            }
            return solenoid;

        }
        else {
            free_solenoid(solenoid);
            return NULL;
        }

    }
    else {
        printf("Error: failed to allocate memory for solenoid.\n");
        return NULL;
    }
}

// Free memory for a Solenoid
void free_solenoid(Solenoid *solenoid) {
    if (solenoid != NULL) {
        free(solenoid->x); free(solenoid->y); free(solenoid->z); 
        free(solenoid->vol); 
        free(solenoid->Jx); free(solenoid->Jy); free(solenoid->Jz); 
        free(solenoid);
    }
}



// --- 
// Outputs
// ---


Line *new_line(Axis axis, double start, double end, uint32_t n) {
    Line *line = calloc(1, sizeof(Line));
    if (line != NULL) {
        line->n = n; 

        line->x = zeros(n); line->y = zeros(n); line->z = zeros(n); 
        line->Bx = zeros(n); line->By = zeros(n); line->Bz = zeros(n); 

        if ( (line->x != NULL) && (line->y != NULL) && (line->z != NULL) && \
            (line->Bx != NULL) && (line->By != NULL) && (line->Bz != NULL)) {
                
            // TODO: this does not error if end < start, is that ok?
            double step = (end - start) / (double)(n - 1);

            switch (axis) {
                case X: 
                    line->x[0] = start;
                    for (uint32_t i=1; i<n; i++) {
                        line->x[i] = line->x[i-1] + step;
                    }
                    break;

                case Y: 
                    line->y[0] = start;
                    for (uint32_t i=1; i<n; i++) {
                        line->y[i] = line->y[i-1] + step;
                    }
                    break;

                case Z: 
                    line->z[0] = start;
                    for (uint32_t i=1; i<n; i++) {
                        line->z[i] = line->z[i-1] + step;
                    }
                    break;
                }
            return line; 
        } 
        else {
            free_line(line);
            return NULL;
        }

    }
    else {
        printf("Error in allocating memory for the line().\n");
        return NULL;
    }
}


void free_line(Line *line) {
    if (line != NULL) {
        free(line->x);
        free(line->y);
        free(line->z);
        free(line->Bx);
        free(line->By);
        free(line->Bz);
    }
    free(line);
}