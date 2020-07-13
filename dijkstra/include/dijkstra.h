#ifndef __DIJKSTRA__
#define __DIJKSTRA__

#define infty 9999
#include "binheap.h"
#include <stdio.h>

typedef struct queue{  // array based queue structure
    int* array;
    size_t size;
}queue;

// node in adjacency list
typedef struct adjnode{
    int weight;       // edge weight
    int id;           // identifier for the node
}adjnode;

typedef struct graph{
    size_t dim;     // number of nodes in the graph
    int* V;         // array in which are randomly inserted the number of neighbours of the node whose id corresponds to the position in V of a specific number of neighbours
    adjnode** adj;  // graph as an array of adjacency lists
}graph;

typedef struct pair{ // type returned by dijkstra in order to verify the correctness of the algorithm
    int* dist;       // array of distances between the source and all nodes
    int* pred;       // array of predecessors
}pair;



graph* build_graph(size_t dim);

void deallocate_graph(graph* G);

int isempty_array(queue* q);

int extract_min_array(queue* q, int* dist);

int extract_min_heap(binheap_type* q);

void init(graph *G, int* dist, int* pred);

queue* build_queue(size_t dim);

void swap(int* arr, int i, int j);

void relax_distance_array(int i, int j, int weight, int* dist, int* pred);

void relax_distance_heap(binheap_type* q, int i, int j, int weight, int* dist, int* pred);

void update_distance_array(int i, int new_dist, int* dist);

void update_distance_heap(binheap_type* q, int i, const void* new_dist);

pair dijkstra(graph *G, int source);

pair dijkstra_heap(graph *G, int source);

#endif
