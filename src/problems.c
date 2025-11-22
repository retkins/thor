#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

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

        loop->x = NULL; loop->y = NULL; loop->z = NULL; loop->vol = NULL; 
        loop->Jx = NULL; loop->Jy = NULL; loop->Jz = NULL;

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
            // communication of errors is done in `zeros()`
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



// --- 
// Outputs
// ---

// A line in 3D space drawn along the specified axis
// typedef struct Line {
//     uint32_t n;                 // number of points
//     double *x, *y, *z;          // [m]
//     double *Bx, *By, *Bz        // [T]
// } Line;

Line *new_line(Direction dir, double start, double end, uint32_t n) {
    Line *line = calloc(1, sizeof(Line));
    if (line != NULL) {
        line->n = n; 
        
        // Start with null pointers in case we have an alloc failure and need 
        // to free them all at once (for conciseness)
        line->x = NULL; line->y = NULL; line->z = NULL; 
        line->Bx = NULL; line->By = NULL; line->Bz = NULL; 

        line->x = zeros(n); line->y = zeros(n); line->z = zeros(n); 
        line->Bx = zeros(n); line->By = zeros(n); line->Bz = zeros(n); 

        if ( (line->x != NULL) && (line->y != NULL) && (line->z != NULL) && \
            (line->Bx != NULL) && (line->By != NULL) && (line->Bz != NULL)) {
                
            // TODO: this does not error if end < start, is that ok?
            double step = (end - start) / (double)(n - 1);

            switch (dir) {
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
