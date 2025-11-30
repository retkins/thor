#include "../include//internal/sort.h"
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

typedef uint32_t u32; 
typedef uint64_t u64; 

// Naive bubble sort implementation has O(N^2) time complexity
uint64_t *bubble_sort(uint64_t *array, uint32_t n) {
    uint64_t *indices = malloc(n*sizeof(uint64_t));
    for (uint32_t i=0; i<n; i++) {
        indices[i] = i;
    }

    uint64_t temp_value; 
    uint32_t temp_index;
    uint64_t iterations = 0; 
    uint64_t swaps = 0; 
    while (true) {
        bool did_swap = false; 
        for (uint32_t i=0; i<(n-1); i++) {
            if (array[i+1] < array[i]) {
                did_swap = true; 
                // Swap the values first 
                temp_value = array[i]; 
                array[i] = array[i+1];
                array[i+1] = temp_value; 

                // Now swap the indices 
                temp_index = indices[i]; 
                indices[i] = indices[i+1]; 
                indices[i+1] = temp_index;
                swaps++;
            }
        }
        iterations++;
        if (!did_swap) {
            break;
        } 
    }
    printf("Sorted %i-length list using %lli swaps in %lli iterations.\n", n, swaps, iterations);

    return indices;
}


// Comparator function for qsort()
// from https://www.geeksforgeeks.org/c/qsort-function-in-c/#
int comp_int(const void *a, const void *b) {
    return *(int *)b - *(int *)a;     // yeesh!
}

typedef struct IndexValuePair {
    uint32_t index; 
    uint64_t value; 
} IndexValuePair; 

// Comparator 
int comp_pair(const void *a, const void *b) {
    return (*(IndexValuePair *)a).value - (*(IndexValuePair *)b).value;     // omg
}


// TODO: update this with an efficient version
u64 *quick_sort(u64 *array, u64 n) {

    // Setup arrays containing the data to sort 
    IndexValuePair *pairs = malloc(n*sizeof(IndexValuePair)); 
    u64 *indices = malloc(n*sizeof(u64));
    for (u32 i=0; i<n; i++) {
        pairs[i].index = i; 
        pairs[i].value = array[i];
    }

    // In a sane world, we'd do this:
    // qsort(array, n, sizeof(uint64_t),comp_int);
    // but we need the indexes
    qsort(pairs, n, sizeof(IndexValuePair), comp_pair);

    // Bring the values back into the two arrays 
    for (u32 i=0; i<n; i++) {
        array[i] = pairs[i].value; 
        indices[i] = pairs[i].index;
    }

    free(pairs);
    return indices;
}

// TODO: implement Radix Sort 
u64 *radix_sort(u64 *array, u64 n) {
    u64 *indices;
    return indices;
}

uint64_t *sort(uint64_t *array, uint32_t n, SortMethod method) {

    if (method == BubbleSort) {
        return bubble_sort(array, n);
    }
    else if (method == QuickSort) {
        return quick_sort(array, n); 
    }
    else {
        // return radix_sort(array, n);
        return quick_sort(array, n);
    }
}

