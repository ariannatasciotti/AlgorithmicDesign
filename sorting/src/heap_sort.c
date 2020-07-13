#include "heap_sort.h"
#include <binheap.h>
#include <string.h>

void heap_sort(void *A, const unsigned int n, const size_t elem_size, total_order leq){

    binheap_type* h = build_heap(A, n, n, elem_size, leq);
    
    for(unsigned int i = n-1; i > 0; i--) extract_min(h); // from the highest to the lowest 
    
    for(size_t i = 0; i < n/2; i++) swap(A+(i*elem_size), A+((n-i-1)*elem_size), elem_size); // from the lowest to the highest
        
    free(h);
}
