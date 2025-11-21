/*  Setup dicrete problems useful for running tests or evaluating standard 
    electromagnetic conditions.
*/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "../include/internal/utils.h"
#include "../include/internal/problems.h"


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

