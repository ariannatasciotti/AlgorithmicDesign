#include "selection_sort.h"

void selection_sort(void *A, const unsigned int n, const size_t elem_size, total_order leq)
{
    for(size_t i = n-1; i > 0; i--){
        size_t max = 0;
           for(size_t j = 1; j < i+1; j++){
               if(!leq(A+(j*elem_size), A+(max*elem_size))){
                   max = j;
               }
           }
        swap(A+(i*elem_size), A+(max*elem_size), elem_size);
    }
}



