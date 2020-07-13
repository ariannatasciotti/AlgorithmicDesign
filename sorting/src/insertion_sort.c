#include "insertion_sort.h"


void insertion_sort(void *A, const unsigned int n, const size_t elem_size, total_order leq)
{
    size_t i;

    for(size_t j = 1; j < n; j++){
        i = j;
        while (i > 0 && leq(A+(i*elem_size), A+(elem_size*(i-1)))){
            swap(A+(elem_size*(i-1)), A+(elem_size*i), elem_size);
            i--;
        }
    }
}

