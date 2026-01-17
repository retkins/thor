#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h> 

#include "../include/internal/utils.h"
#include "../include/internal/biotsavart.h"
#include "../include/internal/sources.h"


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
// Dense Solenoid
// ---

// typedef struct DenseSolenoid {
//     double inner_radius;            // [m]
//     double outer_radius;            // [m]
//     double length;                  // [m]
//     double zcentroid;               // [m]
//     double current;                 // [A]
//     double turns_per_unit_length;   // [turns/m]
//     uint32_t n;                 // number of points
//     double *x, *y, *z;          // [m] location of source points in loop 
//     double *vol;                // [m^3] volume of each source point 
//     double *Jx, *Jy, *Jz;       // [A/m^2] current density vector at each source point

// } DenseSolenoid;

// Allocate memory for a new DenseSolenoid
DenseSolenoid *new_dense_solenoid(
    double inner_radius, double outer_radius, double length, double zcentroid, 
    double current, double element_size
) {
    DenseSolenoid *ds = calloc(1, sizeof(DenseSolenoid));

    if (ds != NULL) {

        // Determine how many layers and how many turns per layer there area
        int nlayers = (int)ceil(length / element_size); 
        int nturns_per_layer = (int)ceil((outer_radius - inner_radius) / element_size);
        
        // Now step through each turn in a layer to count the number of elements
        int n_elements_per_layer = 0;
        int *n_elements_per_radius = malloc(nturns_per_layer*sizeof(int));
        double r = inner_radius + 0.5*element_size; 
        double rstep = (outer_radius - inner_radius) / (double)nturns_per_layer;

        for (int i=0; i<nturns_per_layer; i++) {
            int n_elements = (int)ceil(2*PI*r/element_size); 
            n_elements_per_layer += n_elements; 
            n_elements_per_radius[i] = n_elements;
            r += rstep;
        }

        // That was pretty gnarly... 

        // Total number of elements 
        int n = nlayers * n_elements_per_layer;
        int nturns = nturns_per_layer * nlayers;

        ds->inner_radius = inner_radius; ds->outer_radius = outer_radius;
        ds->length = length; ds->zcentroid = zcentroid; ds->current = current; 
        ds->element_size = element_size; ds->n = n;
        ds->turns_per_unit_length = (double)nturns / length;

        ds->x = zeros(n); ds->y = zeros(n); ds->z = zeros(n); ds->vol = zeros(n);
        ds->Jx = zeros(n); ds->Jy = zeros(n); ds->Jz = zeros(n); 

        if ( (ds->x != NULL) && (ds->y != NULL) && (ds->z != NULL) && \
            (ds->vol != NULL) && (ds->Jx != NULL) && (ds->Jy != NULL) && 
            (ds->Jz != NULL) 
        ) {

            double element_volume = pow(element_size, 3.0); 
            double element_area = element_size*element_size;
            double element_jdensity = current/element_area;

            double z = zcentroid - length / 2.0; 
            double zstep = length / (double)(nlayers - 1); 

            int e = 0;  // counter for all elements 

            for (int i=0; i<nlayers; i++) {
                // outer loop over layers 

                double r = inner_radius; 
                for (int j=0; j<nturns_per_layer; j++) {
                    // middle loop over consecutive radii 
                    // int n_per_radius = (int)ceil(2*PI*r / element_size);
                    double theta_step = 2*PI/(double)n_elements_per_radius[j];
                    double theta = 0.0; 
                    for (int k=0; k<n_elements_per_radius[j]; k++) {
                        // inner loop over elements in a radius
                        ds->x[e] = r * cos(theta); 
                        ds->y[e] = r * sin(theta); 
                        ds->z[e] = z;
                        ds->vol[e] = element_volume;
                        ds->Jx[e] = -element_jdensity * sin(theta); 
                        ds->Jy[e] = element_jdensity *cos(theta); 
                        e++;
                        if (e > n) {
                            printf("Error; e > n\n");
                        }
                        theta += theta_step;
                    }
                    r += element_size;
                }
                z += zstep;
            }
            return ds;

        }
        else {
            free_dense_solenoid(ds);
            return NULL;
        }
    }
    else {
        fprintf(stderr, "Error allocating memory for Dense Solenoid.\n");
        return NULL;
    }
}

// Deallocate memory for a DenseSolenoid
void free_dense_solenoid(DenseSolenoid *ds) {
    if (ds != NULL) {
        free(ds->x); free(ds->y); free(ds->z); free(ds->vol);
        free(ds->Jx); free(ds->Jy); free(ds->Jz);
        free(ds);
    }
}