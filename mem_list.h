#ifndef MEM_LIST_H
#define MEM_LIST_H
#include <inttypes.h>
#include "vector.h"

typedef struct tri_mem_list
{
  unsigned char  *** t_arr; //Triangle array
  unsigned short ** sym_index; //Convert vertex to the index of it's symmetry class. First dimension is vertex_to_index, second is symmetrys
  arr3 dim_size; //Size of the 3D axis
  arr2 dim_mult; //Multiplies for the vertex_to_index
  arr3 dim;      //Original dimension
  int  fund;     //if true, we only store the triangles with a vertex in the fund domain
  
} tri_mem_list;

typedef unsigned short tri_index[3];

typedef struct tri_index_list
{
  tri_index * index_list;
  size_t len;
  arr3 dim;
} tri_index_list;



#define vertex_to_index(vertex,dim_mult) (vertex[0] * dim_mult[0] + vertex[1] * dim_mult[1] + vertex[2])
//Set mem_list index
#define SMI(t_arr,index) (t_arr[index[0]][index[1]][index[2] / 8] |= 1 << (index[2] % 8))
//Returns whether the index bit is set
#define GMI(t_arr,index) (t_arr[index[0]][index[1]][index[2] / 8] & (1 << (index[2] % 8)))
//Exist mem_list index
#define EMI(t_arr,index) ((t_arr[index[0]][index[1]]))
//Clear mem_list index
#define CMI(t_arr,index) (t_arr[index[0]][index[1]][index[2] / 8] &= ~(1 << (index[2] % 8)))

#define vertices_to_index_fund(v1,v2,v3,mem_list,indices) (\
          vertices_unique_fund(vertex_to_index(v1,mem_list->dim_mult),\
                               vertex_to_index(v2,mem_list->dim_mult),\
                               vertex_to_index(v3,mem_list->dim_mult),mem_list,indices)) 


int index_list_from_file(tri_index_list * result, char * filename);
int index_list_to_file(tri_index_list * list, char * filename);


tri_mem_list mem_list_from_triangle_list(triangle_list * list);
tri_mem_list mem_list_from_index_list(tri_index_list * list);
tri_mem_list mem_list_from_index_list_fund(tri_index_list * list);

int mem_list_from_file(tri_mem_list * result, char * filename);
int mem_list_to_file(tri_mem_list * list, char * filename);
triangle_list mem_list_to_triangle_list(tri_mem_list * list);

size_t mem_list_count(tri_mem_list * list);
size_t mem_list_memory(tri_mem_list * list);
size_t mem_list_indices(tri_mem_list * list, tri_index_list * index_list);
tri_mem_list mem_list_init(arr3 dim, arr3 dim_size);
tri_mem_list mem_list_init_fund(int dim);
void mem_list_free(tri_mem_list * list);
void mem_list_clean(tri_mem_list * list);


void mem_list_clear(tri_mem_list * list, ptriangle triang);
void mem_list_set(tri_mem_list * list, ptriangle triang);
int  mem_list_get(tri_mem_list * list, ptriangle triang);

void mem_list_clear_sym(tri_mem_list * list, ptriangle triang);
void mem_list_set_sym(tri_mem_list * list, ptriangle triang);

void mem_list_clear_sym_fund(tri_mem_list * list, ptriangle triang);
void mem_list_set_sym_fund(tri_mem_list * list, ptriangle triang);
int mem_list_get_sym_fund(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3);


int index_in_array(int index, int * indices, size_t len);
void vertices_to_index(arr3 v1, arr3 v2, arr3 v3, arr2 dim_mult, tri_index indices);
void triangle_to_index(ptriangle triang, arr2 dim_mult, tri_index indices);
int vertices_unique_fund(unsigned short idx1, unsigned short idx2, unsigned short idx3, tri_mem_list * mem_list, tri_index indices);

int triangle_to_index_fund(ptriangle triang, tri_mem_list * mem_list, tri_index indices);

triangle triangle_from_index(tri_index indices, arr3 dim);
void vertex_from_index(unsigned short index, arr2 dim_mult,arr3 vertex);

#endif
