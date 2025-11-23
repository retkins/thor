#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h> 

#include "../include/internal/utils.h"
#include "../include/internal/biotsavart.h"
#include "../include/internal/targets.h"

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