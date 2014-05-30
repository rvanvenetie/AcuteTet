#ifndef MEM_LIST_H
#define MEM_LIST_H
#include <inttypes.h>
#include "vector.h"


typedef struct tri_mem_tet{
  vert_index_array vert_to_index; 
  int tet_len; //Amount of vertices in the tetrahedron
} tri_mem_tet;

typedef struct tri_mem_cube {
  arr3 dim_size; //Size of the data 
  arr2 dim_mult;
} tri_mem_cube;

typedef struct tri_mem_fund {
  int dim_size; //Size of the data 
  arr2 dim_mult; //Multiplies for the vertex_to_index
  vert_index ** sym_index; //Cached symmetries on each vertex index. First dimension is the index, second the symmetry number
} tri_mem_fund;

typedef struct tri_mem_fund2 {
  vert_index_array vert_to_index; //Convert point to it's index
  arr3 *           vert_from_index; //Convert index to it's vertex
  vert_index **  sym_index;     //Cached symmetries
  size_t cube_len; //Total amount of points in the cube
  size_t fund_len; //Amount of points inside the fundamental domain
} tri_mem_fund2;

typedef struct tri_mem_list
{
  unsigned char  *** t_arr; //Actual data array
  union {
    tri_mem_cube mem_cube;
    tri_mem_fund mem_fund;
    tri_mem_tet mem_tet;
    tri_mem_fund2 mem_fund2;
  };
  arr3 dim;      //Original dimension
  int  mode;     //Stores the type of mem_list (store every triangle in the cube, only in fundamental domain, or triangles in the unit tetrahedron)
  
} tri_mem_list;

#define MEM_LIST_CUBE  0
#define MEM_LIST_FUND  1
#define MEM_LIST_FUND2 2
#define MEM_LIST_TET   3
typedef unsigned short tri_index[3];
typedef unsigned short tet_index[4];

typedef struct tri_index_list
{
  tri_index * index_list;
  size_t len;
  arr3 dim;
} tri_index_list;


//Set mem_list index
#define SMI(t_arr,index) (t_arr[index[0]][index[1]][index[2] / 8] |= 1 << (index[2] % 8))
//Returns whether the index bit is set
#define GMI(t_arr,index) (t_arr[index[0]][index[1]][index[2] / 8] & (1 << (index[2] % 8)))
//Exist mem_list index
#define EMI(t_arr,index) ((t_arr[index[0]][index[1]]))
//Clear mem_list index
#define CMI(t_arr,index) (t_arr[index[0]][index[1]][index[2] / 8] &= ~(1 << (index[2] % 8)))


/*
 * Macros to convert a vertex to it's index in the mem_list
 */
#define vertex_to_index_tet(vertex, vertex_index_array) (vertex_index_array[vertex[0]][vertex[1]][vertex[2]])
#define vertex_to_index_cube(vertex,dim_mult) (vertex[0] * dim_mult[0] + vertex[1] * dim_mult[1] + vertex[2])
#define vertex_to_index_fund(vertex,dim_mult) (vertex[0] * dim_mult[0] + vertex[1] * dim_mult[1] + vertex[2])
#define vertex_to_index_fund2(vertex, vertex_index_array) (vertex_index_array[vertex[0]][vertex[1]][vertex[2]])

/* 
 * See indices_unique_*. This is a wrapper to convert vertices to a unique index directly
 */
#define vertices_to_index_fund(v1,v2,v3,mem_list,indices) (\
          indices_unique_fund(vertex_to_index_cube(v1,mem_list->mem_fund.dim_mult),\
                              vertex_to_index_cube(v2,mem_list->mem_fund.dim_mult),\
                              vertex_to_index_cube(v3,mem_list->mem_fund.dim_mult),mem_list,indices)) 
                               
#define vertices_to_index_cube(v1, v2, v3, dim_mult, indices) (\
          indices_unique_cube(vertex_to_index_cube(v1, dim_mult),\
                              vertex_to_index_cube(v2, dim_mult),\
                              vertex_to_index_cube(v3, dim_mult),indices))

#define vertices_to_index_tet(v1, v2, v3, vert_to_index, indices) (\
          indices_unique_tet(vertex_to_index_tet(v1, vert_to_index),\
                             vertex_to_index_tet(v2, vert_to_index),\
                             vertex_to_index_tet(v3, vert_to_index),indices))
                             
/*
 * These are wrappers to convert triangles (thus three vertices) to a unique index directly
 */
#define triangle_to_index_fund(triang, mem_list, indices) (\
          vertices_to_index_fund(triang.vertices[0],triang.vertices[1], triang.vertices[2], mem_list,indices))
          
#define triangle_to_index_cube(triang, dim_mult, indices) (\
          vertices_to_index_cube(triang.vertices[0],triang.vertices[1], triang.vertices[2], dim_mult,indices))
          
#define triangle_to_index_tet(triang, vert_to_index, indices) (\
          vertices_to_index_tet( triang.vertices[0],triang.vertices[1], triang.vertices[2], vert_to_index,indices))



/*
 * Index - Vertex functions
 */
void indices_unique_cube(vert_index idx1, vert_index idx2, vert_index idx3, tri_index indices);
void indices_unique_tet(vert_index idx1, vert_index idx2, vert_index idx3, tri_index indices);
int indices_unique_fund(vert_index idx1, vert_index idx2, vert_index idx3, tri_mem_list * mem_list, tri_index indices);

/*
 * Functions that generate vertex -> index arrays
 */
void gen_vertex_to_index_fund(int dim, vert_index_array * vertex_index, size_t * fund_len, size_t * cube_len);
void gen_vertex_from_index_fund(int dim, arr3 ** index_vertex, vert_index_array vertex_index);
void gen_vertex_to_index_tet(int dim, vert_index_array * vertex_index, int *tet_len);


void vertex_from_index_cube(vert_index index, arr2 dim_mult,arr3 vertex) ;
triangle triangle_from_index_cube(tri_index indices, arr2 dim_mult);
triangle triangle_from_index_fund2(tri_index indices,arr3 * index_vertex);


/*
 * mem_list init/free functions
 */
#define MEM_LIST_FALSE 0
#define MEM_LIST_TRUE 1
tri_mem_list mem_list_init_cube(arr3 dim);
tri_mem_list mem_list_init_fund(int dim, int init_value);
tri_mem_list mem_list_init_fund2(int dim, int init_value);
tri_mem_list mem_list_init_tet(int dim, int init_value);
void mem_list_free(tri_mem_list * list);



/*
 * mem_list operations
 */
int mem_list_dim_size(tri_mem_list * list, int dim, int idx1, int idx2);
int mem_list_row_empty(tri_mem_list * list,int i, int j);
size_t mem_list_count(tri_mem_list * list);
size_t mem_list_memory(tri_mem_list * list);
void mem_list_clean(tri_mem_list * list);
size_t mem_list_indices(tri_mem_list * list, tri_index_list * index_list);
triangle_list mem_list_to_triangle_list(tri_mem_list * list);


/*
 * mem_list file functions
 */
#define MEM_LIST_SAVE_FULL 0
#define MEM_LIST_SAVE_CLEAN 1
int mem_list_from_file(tri_mem_list * result, char * filename);
int mem_list_to_file(tri_mem_list * list, char * filename, int mode);









/*
 * mem_list set/get/clear cube functions
 */
void mem_list_clear_cube(tri_mem_list * list, ptriangle triang);
void mem_list_set_cube(tri_mem_list * list, ptriangle triang);
int  mem_list_get_cube(tri_mem_list * list, ptriangle triang);
void mem_list_clear_sym_cube(tri_mem_list * list, ptriangle triang);
void mem_list_set_sym_cube(tri_mem_list * list, ptriangle triang);

/*
 * mem_list set/get/clear fund functions
 */
void mem_list_clear_fund(tri_mem_list * list, ptriangle triang);
void mem_list_set_fund(tri_mem_list * list, ptriangle triang);
int mem_list_get_fund(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3);
/*
 * mem_list set/get/clear tet functions
 */
void mem_list_clear_fund2(tri_mem_list * list, ptriangle triang);
void mem_list_set_fund2(tri_mem_list * list, ptriangle triang);
int mem_list_get_fund2(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3);

/*
 * mem_list set/get/clear tet functions
 */
int mem_list_get_tet(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3);
void mem_list_set_tet(tri_mem_list * list, ptriangle triang);
void mem_list_clear_tet(tri_mem_list * list, ptriangle triang);
/*
 * old
 */
tri_mem_list mem_list_from_index_list(tri_index_list * list);
tri_mem_list mem_list_from_index_list_fund(tri_index_list * list);
int index_list_from_file(tri_index_list * result, char * filename);
int index_list_to_file(tri_index_list * list, char * filename);


#endif
