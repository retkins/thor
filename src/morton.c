
#include "../include/internal/morton.h"
#include "../include/internal/sort.h"

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>


BBox new_bbox(double xc, double yc, double zc, double side_length) {
    double half = side_length/2.0;
    return (BBox){
        .xc=xc, .yc=yc, .zc=zc, .side_length=side_length,
        .xmin=xc-half, .ymin=yc-half, .zmin=zc-half
    };
}


// Determine the scale factor based on the maximum quantized depth of the tree
double calculate_scale_factor(uint32_t L) {
    return (double)(1u << L);       // TODO: should this be `1u << L` or `(1u <<L) - 1`?
}


// Normalize a point's location in a grid to the range [0,1] (inclusive)
void normalize(double normalized_point[3], double point[3], BBox bbox) {
    normalized_point[0] = (point[0] - bbox.xmin)/bbox.side_length;
    normalized_point[1] = (point[1] - bbox.ymin)/bbox.side_length;
    normalized_point[2] = (point[2] - bbox.zmin)/bbox.side_length;
}


// Quantize a point's location into an integer 
void quantize(uint32_t quantized_point[3], double normalized_point[3], double scale) {
    quantized_point[0] = (uint32_t)floor(scale * normalized_point[0]);
    quantized_point[1] = (uint32_t)floor(scale * normalized_point[1]);
    quantized_point[2] = (uint32_t)floor(scale * normalized_point[2]);
}


// Interleave the bits of three integers into a Morton code 
// Reference: https://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/
uint64_t interleave(uint32_t quantized_point[3]) {
    uint64_t code = 0;
    uint32_t x = quantized_point[0]; 
    uint32_t y = quantized_point[1];
    uint32_t z = quantized_point[2];

    for (uint64_t i=0; i<21; i++) {
        code |= ((x & ((uint64_t)1 << i)) << 2*i) | ((y & ((uint64_t)1 << i)) << (2*i + 1)) | ((z & ((uint64_t)1 << i)) << (2*i + 2));
    }
    return code;
}


// Convert a position into a morton code 
uint64_t encode(double point[3], double scale, BBox bbox) {
    double normalized_point[3]; 
    uint32_t quantized_point[3]; 
    normalize(normalized_point, point, bbox); 
    quantize(quantized_point, normalized_point, scale);
    return interleave(quantized_point);
}


uint64_t *encode_points(double *x, double *y, double *z, uint32_t n, uint32_t L, BBox bbox) {
    uint64_t *codes = calloc(n, sizeof(uint64_t)); 
    double scale = calculate_scale_factor(L);

    double point[3]; 
    for (uint32_t i=0; i<n; i++) {
        point[0] = x[i]; 
        point[1] = y[i]; 
        point[2] = z[i]; 
        codes[i] = encode(point,scale, bbox);
    }

    return codes;
}


// --- 
// Tests 
// --- 


// Create a random array of integers
uint64_t *rand_ints(uint32_t n, uint64_t max) {
    srand((uint64_t)clock());
    uint64_t *array = calloc(n, sizeof(uint64_t)); 
    for (uint32_t i=0; i<n; i++) {
        array[i] = max*rand()/RAND_MAX;
    }

    return array;
}

// Print an array of integers to the commande lin
void print_ints(uint64_t *array, uint32_t n) {
    for (uint32_t i=0; i<n; i++) {
        printf("%li ", array[i]);
    }
    printf("\n");
}


void test_interleave() {
    // From the above reference 
    uint32_t point[3] = {5, 9, 1}; 
    uint64_t expected_result = 1095; 

    uint64_t code = interleave(point);
    printf("Testing interleaving function...");
    if (code == expected_result) { printf("Test PASSED\n"); }
    else { printf("Test FAILED.\n"); }
    printf("\tMorton code for (%i, %i, %i) = %li\n", point[0], point[1], point[2], code);
}


void test_encode(uint32_t L) {
    BBox bbox = new_bbox(0.5, 0.5, 0.5, 1.0);
    double point[3] = {0.75, 0.75, 0.75};
    uint64_t code = encode(point, L, bbox); 
    printf("Testing Morton encoding...\n");
    printf("\tCode = %li\n", code);
}


void test_sort(uint32_t n) {
    // 
    uint64_t *array = rand_ints(n, 10*n);
    // printf("Before sorting:\n");
    // print_ints(array, n);
    printf("Bubble sort:\n");
    clock_t start = clock();
    sort(array, n, BubbleSort);
    clock_t end = clock(); 
    printf("Elapsed time: %.6f\n", (double)(end - start)/CLOCKS_PER_SEC);
    free(array); 
    array = rand_ints(n, 10*n);

    printf("qsort:\n");
    start = clock() ;
    sort(array, n, QuickSort);
    end = clock();
    printf("Elapsed time: %.6f\n", (double)(end - start)/CLOCKS_PER_SEC);
    // printf("After sorting:\n"); 
    // print_ints(array, n);
    free(array);
}

// int main() {
//     test_interleave();
//     test_encode(3);

//     // test_sort(1000);
//     return 0;
// }