#include "dijkstra.h"
#include <stdio.h>
#include <stdlib.h>

graph* build_graph(size_t dim)
{
    graph* G = (graph*)malloc(sizeof(graph));
    G->dim = dim;
    G->V = malloc(sizeof(int)*dim);
    G->adj = malloc(sizeof(adjnode*)*dim);
    for (size_t i = 0; i < dim; i++){
        G->V[i] = rand()%5;
        G->adj[i] = malloc(sizeof(adjnode)*G->V[i]);
        for(size_t j = 0; j < G->V[i]; j++){
            G->adj[i][j].id = rand()%dim;     // id for the neighbours of node i
            G->adj[i][j].weight = 1+rand()%10; // weight of the edge between node i and its neighbour j
            if(i == G->adj[i][j].id){         // to avoid self loops
                G->adj[i][j].id = rand()%dim;
            }
            //printf("starting node:%zu, destination node:%d, weight:%d\n", i, G->adj[i][j].id,  G->adj[i][j].weight);
        }
    }
    return G;
}


void deallocate_graph(graph* G)
{
    free(G->V);
    for(size_t i = 0; i < G->dim; i++){
       free(G->adj[i]);
     }
    free(G->adj);
    free(G);
}


int isempty_array(queue* q)
{
    return q->size == 0;
}


int extract_min_array(queue* q, int* dist)
{
    size_t index = 0;
 
    for(size_t i = 1; i < q->size; i++){
        if(dist[q->array[i]] < dist[q->array[index]]){
            index = i;
        }
    }
    swap(q->array, q->size-1, index);
    q->size--;
    //printf("%d\n", q->size);
    return q->size;
}


int extract_min_heap(binheap_type* q)
{
    const void* n = extract_min(q);
    // index in A, need an int for my implementation
    return (n - q->A)/q->key_size;
}


void relax_distance_array(int id_u, int id_v, int weight, int* dist, int* pred)
{
    if(dist[id_u] + weight < dist[id_v]){
        update_distance_array(id_v, dist[id_u]+weight, dist);
        pred[id_v] = id_u;
    }
}


void update_distance_array(int id_u, int new_dist, int* dist)
{
    dist[id_u] = new_dist;
}


void relax_distance_heap(binheap_type* q, int id_u, int id_v, int weight, int* dist, int* pred)
{
    int* dist_temp = q->A;    // heap is constructed on the array of distances
    if((dist_temp[id_u] + weight) < dist_temp[id_v]){
        int temp = dist_temp[id_u] + weight;
        update_distance_heap(q, id_v, (void*) &temp);
        pred[id_v] = id_u;
    }
}


void update_distance_heap(binheap_type* q, int id_u, const void* new_dist)
{
    void * ptr_id = q->A + (id_u)*(q->key_size);
    decrease_key(q, ptr_id, new_dist);
}


void init(graph *G, int* dist, int* pred)
{
    for(size_t i = 0; i < G->dim; i++){
        dist[i] = infty;    // the arrays of the distances and the predecessors are allocated inside the dijkstra function
        pred[i] = -1;
    }
}

queue *build_queue(size_t dim)
{
    queue *queue = malloc(sizeof(queue));
    queue->size = dim;
    queue->array = malloc(sizeof(int)*dim);
    for(size_t i = 0; i < dim; i++)
        queue->array[i] = i;
    return queue;
}


pair dijkstra(graph* G, int source)
{
    int* dist = malloc(sizeof(int)*G->dim);
    int* pred = malloc(sizeof(int)*G->dim);
    
    init(G, dist, pred);
    
    dist[source] = 0;  // distance of source from source is 0
    
    queue* q = build_queue(G->dim);
    
    while(!isempty_array(q))
    {
        int v = extract_min_array(q, dist);
        //printf("%d\n", q->array[v]);
        
        for(size_t i = 0; i < G->V[q->array[v]]; i++){
            //printf("%d %d %d\n", q->array[v], G->adj[q->array[v]][i].id, G->adj[q->array[v]][i].weight);
            relax_distance_array(q->array[v], G->adj[q->array[v]][i].id,  G->adj[q->array[v]][i].weight, dist, pred);
            //printf("%d %d %d\n", q->array[v], G->adj[q->array[v]][i].id, G->adj[q->array[v]][i].weight);
        }
        
    }
    
    free(q->array);
    free(q);
    pair results;
    results.dist = dist;
    results.pred = pred;
    return results;
}


void swap(int* arr, int i, int j)
{
    int temp;
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}

   
pair dijkstra_heap(graph* G, int source)
{
    int* dist = malloc(sizeof(int)*G->dim);
    int* pred = malloc(sizeof(int)*G->dim);
    
    init(G, dist, pred);
    
    dist[source] = 0;
    
    binheap_type* q = build_heap(dist, G->dim, G->dim, sizeof(int), leq_int);    // build the heap on the array of the distances
    
    while(!is_heap_empty(q))
    {
        int v = extract_min_heap(q);
        //printf("%d\n", v);
        
        for(size_t i = 0; i < G->V[v]; i++){
            //printf("%d %d %d\n", v, G->adj[v][i].id, G->adj[v][i].weight);
            relax_distance_heap(q, v, G->adj[v][i].id, G->adj[v][i].weight, dist, pred);
            //printf("%d %d %d\n", v, G->adj[v][i].id, G->adj[v][i].weight);
        }
        
    }

    delete_heap(q);
    pair results;
    results.dist = dist;
    results.pred = pred;
    return results;
}
