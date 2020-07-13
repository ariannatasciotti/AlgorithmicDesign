#include <binheap.h>
#include <string.h>
#include <stdio.h> 

#include <binheap.h>
#include <string.h>
#include <stdio.h>

#define PARENT(node) ((node-1)/2)

#define LEFT_CHILD(node) (2*(node)+1)

#define RIGHT_CHILD(node) (2*(node+1))

#define VALID_NODE(H,node) ((H)->num_of_elem>(node))

#define ADDR(H, node) ((H)->A+(node)*((H)->key_size))

#define INDEX_OF(H,addr) (((addr)-((H)->A))/(H)->key_size)

#define KEY(node) (H->key_pos[node]) // position of the key of a node

#define REV(node) (H->rev_pos[node]) // node corresponding to a given position in A


int is_heap_empty(const binheap_type *H)
{
    return H->num_of_elem == 0;
}


const void *min_value(const binheap_type *H)
{
    if(is_heap_empty(H)) return NULL;
   
    // the minimum is stored in the root aka pos_key[0]
    return ADDR(H, KEY(0));
}


void swap_keys(binheap_type *H, unsigned int n_a, unsigned int n_b)
{
    unsigned int tmp = KEY(n_a);
    KEY(n_a) = KEY(n_b);
    KEY(n_b) = tmp;
    
    // also update rev_pos
    tmp = REV(KEY(n_a));
    REV(KEY(n_a)) = REV(KEY(n_b));
    REV(KEY(n_b)) = tmp;
}


void heapify(binheap_type *H, unsigned int node){

    unsigned int dst_node = node, child;

    do
    {
        dst_node = node;
        
        // find the minimum among node and its children
        child = RIGHT_CHILD(REV(node));

        if(VALID_NODE(H,child) && H->leq(ADDR(H, KEY(child)), ADDR(H, dst_node))) dst_node = KEY(child);
                
        child = LEFT_CHILD(REV(node));

        if(VALID_NODE(H,child) && H->leq(ADDR(H, KEY(child)), ADDR(H, dst_node))) dst_node = KEY(child);
        
        // if the minimum is not in node, swap the keys
        if(dst_node != node) swap_keys(H, REV(dst_node), REV(node));
        
    }while(dst_node != node);
}


const void *extract_min(binheap_type *H)
{
    if(is_heap_empty(H)) return NULL;

    // swapping the keys among the root and the right most leaf of the last level (A[num_of_elem-1])
    swap_keys(H, 0, H->num_of_elem-1);
    
    // deleting the right most leaf of the last level (A[num_of_elem-1])
    H->num_of_elem--;

    heapify(H,KEY(0));
    
    return ADDR(H, KEY(H->num_of_elem));
}


const void *find_the_max(void *A,const unsigned int num_of_elem,const size_t key_size,total_order_type leq)
{
    if(num_of_elem == 0) return NULL;
    
    const void *max_value = A;
    
    // for all the values in A
    for(const void *addr = A+key_size; addr != A+num_of_elem*key_size; addr += key_size)
    {
        // if addr>max_value
        if(!leq(addr,max_value)) max_value = addr;
    }
    
    return max_value;
}


binheap_type *build_heap(void *A,const unsigned int num_of_elem,const unsigned int max_size,const size_t key_size,total_order_type leq)
{
    binheap_type *H = (binheap_type*)malloc(sizeof(binheap_type));
    
    H->A = A;
    H->num_of_elem = num_of_elem;
    H->max_size = max_size;
    H->key_size = key_size;
    H->leq = leq;
    H->max_order_value = malloc(key_size);
    H->key_pos = (int*)malloc(sizeof(int)*max_size);
    H->rev_pos = (int*)malloc(sizeof(int)*max_size);

    
    for(unsigned int i=0; i < max_size; i++){
        H->key_pos[i] = i;
        H->rev_pos[i] = i;
    }
    
    if(num_of_elem == 0) return H;
    
    // get the maximum among A[:num_of_elem-1] and store it in max_order_value
    const void *value = find_the_max(A, num_of_elem, key_size, leq);
    
    memcpy(H->max_order_value, value, key_size);
    
    // fix the heap property from the second last level up to the root
    for(unsigned int i = num_of_elem/2; i > 0; i--) heapify(H, KEY(i));
    
    heapify(H, KEY(0));
    
    return H;
}


void delete_heap(binheap_type *H)
{
    free(H->max_order_value);
    free(H->key_pos);
    free(H->rev_pos);
    free(H);
}


const void *decrease_key(binheap_type *H, void *node, const void *value)
{
    unsigned int node_idx = INDEX_OF(H,node);

    // if node does not belong to H or *value>=*node return null
    if(!VALID_NODE(H, REV(node_idx)) || !(H->leq(value,node))) return NULL;
    
    memcpy(node, value, H->key_size);

    unsigned int parent_idx = KEY(PARENT(REV(node_idx)));
    void *parent = ADDR(H, parent_idx);
    
    while ((REV(node_idx) != 0) && (!H->leq(parent,node)))
    {
        swap_keys(H, REV(parent_idx), REV(node_idx));
        parent_idx = KEY(PARENT(REV(node_idx)));
        parent = ADDR(H,parent_idx);
    }
    
    return node;
}


const void *insert_value(binheap_type *H, const void *value)
{
    // if the heap is already full
    if(H->max_size == H->num_of_elem) return NULL;

    // if the heap is empty or
    if(H->num_of_elem == 0 || !H->leq(value,H->max_order_value)) memcpy(H->max_order_value, value, H->key_size);
    
    // get the position of the new node
    void *new_node_addr = ADDR(H, KEY(H->num_of_elem));
    
    memcpy(new_node_addr, H->max_order_value, H->key_size);
    
    // increase the size of the heap
    H->num_of_elem++;
    
    // decrease the key of the new node
    return decrease_key(H, new_node_addr, value);
}


void print_heap(const binheap_type *H,
                void (*key_printer)(const void *value))
{
    unsigned int next_level_node = 1; // stores the index of the left most node of the next level
    
    for(unsigned int node = 0; node < H->num_of_elem; node++){
        if(node == next_level_node){
            printf("\n");
            next_level_node = LEFT_CHILD(node);
        }else{
            printf("\t");
        }
     key_printer(ADDR(H, KEY(node)));
    }
    printf("\n");
}
