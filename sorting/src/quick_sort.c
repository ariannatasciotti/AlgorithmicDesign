#include "quick_sort.h"
#include "swap.h"

int partition(void *A, size_t i, size_t j, size_t p, const size_t elem_size, total_order leq)
{

    swap(A+(p*elem_size), A+(i*elem_size), elem_size);
    p = i;
    i++;
    j--;

    while(i <= j)
    {
        if(!(leq(A+(i*elem_size), A+(p*elem_size)))){
            swap(A+(i*elem_size), A+(j*elem_size), elem_size);
            j--;
        }
        else{
            i++;
        }
    }
        
    swap(A+(p*elem_size), A+(j*elem_size), elem_size);
    return j;
}


int* partition_3(void *A, size_t i, size_t j, size_t p, const size_t elem_size, total_order leq)
{
    
    int* indexes = malloc(2*sizeof(int));
    swap(A+(p*elem_size), A+(i*elem_size), elem_size);
    p = i;
    i++;
    j--;
    int equal = 0;
    while(i <= j)
    {
        // smaller than pivot
        if(leq(A+(i*elem_size), A+(p*elem_size)) && !((leq(A+(i*elem_size), A+(p*elem_size))) && ((leq(A+(p*elem_size), A+(i*elem_size)))))){
            swap(A+(i*elem_size), A+((p-equal)*elem_size), elem_size);
            p = i;
            i++;
        }
        
         // greater than pivot
         else if(!(leq(A+(i*elem_size), A+(p*elem_size)))){
            swap(A+(i*elem_size), A+(j*elem_size), elem_size);
            j--;
        }
        
        // equal
        else{
            p = i;
            i++;
            equal++;
        }
    }
        
    swap(A+(p*elem_size), A+(j*elem_size), elem_size);
    indexes[0] = j - equal;
    indexes[1] = j;
    return indexes;
}




void quicksort_rec(void *A, size_t l, size_t r, const size_t elem_size, total_order leq)
{
    while(l < r){
        size_t j = partition(A, l, r, l, elem_size, leq);
        quicksort_rec(A, l, j, elem_size, leq);
        l = j+1;
    }
}

void quick_sort(void *A, const unsigned int n,
                const size_t elem_size,
                total_order leq)
{
    quicksort_rec(A, 0, n, elem_size, leq);
}

