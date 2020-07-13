#include "dijkstra.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double test(pair (*f)(graph *, int), graph * G, int source){

    struct timespec requestStart, requestEnd;
    double accum;
    size_t rep = 1;
    
    pair temp_result;
    
    clock_gettime(CLOCK_REALTIME, &requestStart);
    for(size_t i = 0; i < rep; ++i ) temp_result = f(G, source);

    clock_gettime(CLOCK_REALTIME, &requestEnd);

    accum = (requestEnd.tv_sec - requestStart.tv_sec) +
            (requestEnd.tv_nsec - requestStart.tv_nsec) / 1E9;
    
    free(temp_result.dist);
    free(temp_result.pred);
    return accum / rep;
}


int main(void){
    
    size_t n = (1 << 12);
    
    int source = 1;
    
    printf("The source is %d\n", source);
    
    graph* G = build_graph(5);
    
    pair result = dijkstra(G, source);

    printf("Array based version: \n");
    for(size_t i = 0; i < 5; i++){
        printf("node: %zu, distance: %d, predecessor: %d\n", i, result.dist[i], result.pred[i]);
    }
    
    pair result_heap = dijkstra_heap(G, source);
    
    printf("Heap based version: \n");
    for(size_t i = 0; i < 5; i++){
        printf("node: %zu, distance: %d, predecessor: %d\n", i, result_heap.dist[i], result_heap.pred[i]);
    }
    
    free(result.dist);
    free(result.pred);
    free(result_heap.dist);
    free(result_heap.pred);
    deallocate_graph(G);
    
    printf("Performance test: \n");
    printf("Size\tOn Heaps\tOn Arrays\n");
    for (size_t j = 2; j <= n; j *= 2){
        graph* H = build_graph(j);
        printf("%ld\t%f\t%f\n", j, test(dijkstra_heap, H, source),  test(dijkstra, H, source));
        deallocate_graph(H);
    }
    
    return 0;
}
