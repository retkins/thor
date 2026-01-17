/*  Sorting functions 
*/

#ifndef SORT_H 
#define SORT_H

#include <stdint.h>

// Define different sorting algorithms
typedef enum SortMethod {BubbleSort, QuickSort, RadixSort} SortMethod; 

// Sorts an array of 64-bit integers and returns the indices 
uint64_t *sort(uint64_t *array, uint32_t n, SortMethod method);

#endif