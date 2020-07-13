#include "select.h"
#include "quick_sort.h"
#include "swap.h"

unsigned int select_index(void *A, const unsigned int n, const unsigned int i, const size_t elem_size, total_order leq);


unsigned int select_pivot(void* A, const unsigned int n, const size_t elem_size, total_order leq)
{

    if(n < 10){
        quick_sort(A, n, elem_size, leq);
        return n/2;
    }
    
    unsigned int chunks = n/5, c_l, c_r;
    
    for(unsigned int i = 0; i < chunks; i++){
        c_l = i*5 + 1;
        c_r = i*5 + 5;
        quick_sort(A+c_l*elem_size, c_r - c_l, elem_size, leq);
        swap(A+i*elem_size, A+(c_l+2)*elem_size, elem_size);
    }
    
    return select_index(A, chunks-1, chunks/2, elem_size, leq);
}


unsigned int select_index(void *A, const unsigned int n, const unsigned int i, const size_t elem_size, total_order leq)
{
    if(n < 10){
        quick_sort(A, n, elem_size, leq);
        return i;
    }
    
    unsigned int p = select_pivot(A, n, elem_size, leq);
    
    int* k = partition_3(A, 0, n, p, elem_size, leq);
    
    if(i < k[0]) return select_index(A, k[0]-1, i, elem_size, leq);
    
    if(i > k[1]) return select_index(A+k[1]*elem_size, n-k[1]-1, i, elem_size, leq);
    
    return i;
}


void quicksort_aux(void *A, size_t l, size_t r, const size_t elem_size, total_order leq)
{
    while(l < r){
        unsigned int p = l + select_pivot(A+l*elem_size, r-l, elem_size, leq);
        int* k = partition_3(A, l , r, p, elem_size, leq);
        quicksort_aux(A, l, k[0], elem_size, leq);
        l = k[1]+1;
    }
}


void quick_sort_select(void *A, const unsigned int n, const size_t elem_size, total_order leq)
{
   quicksort_aux(A, 0, n, elem_size, leq);
}

