#include <stdlib.h>

typedef struct {
  double value;
  int index;
} pair;


int compare_pairs(const void *pair1, const void *pair2) {
  if ( ((pair*)pair1)->value < ((pair*)pair2)->value )
    return -1;
  else  if ( ((pair*)pair1)->value == ((pair*)pair2)->value )
    return 0;
  else  /* ( ((pair*)pair1)->value > ((pair*)pair2)->value ) */
    return 1;
}

int compare_doubles(const void * a, const void * b)
{
  if ((*(double*)a) == (*(double*)b))
    return 0;
  else if ((*(double*)a) < (*(double*)b))
    return -1;
  else
    return 1;
}


void idxsort_(double* vector, int *sizePtr, int* perm) {
  /*
    vector
    INPUT Array of size double precision elements that are to be sorted in ascending order
    
    size
    INPUT Pointer to the integer number of elements to be sorted.
    An integer in "Fortran side".
    
    perm
    OUTPUT Array of integers, that put vector into ascending order.
    In Matlab notation vector(perm) is an array of values in ascending order.
  */
  
  int i;
  int size = *sizePtr;
  pair* pairs = (pair*) malloc(size * sizeof(pair));
  for(i=0; i < size; i++) {
    pairs[i].value = vector[i];
    pairs[i].index = i;
  }
  qsort(pairs, size, sizeof(pair), compare_pairs);
  
  /* Also take care of the 0- to 1- base convention for the "Fortran side" */
  for(i=0; i < size; i++) {
    perm[i] = pairs[i].index + 1;
  }
  
  /* Finally sort the input vector iself in-place */
  /*  qsort(vector, size, sizeof(int), compare_doubles); */
  qsort(vector, size, sizeof(double), compare_doubles); 
  
}

