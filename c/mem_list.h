#ifndef MEM_LIST_H
#define MEM_LIST_H
#include <inttypes.h>
#include "vector.h"

typedef struct tri_mem_tet{
  vert_index_array vert_to_index; 
  int tet_len; //Amount of vertices in the tetrahedron
} tri_mem_tet;

typedef struct tri_mem_cube {
  int dim_size; //Size of the data 
} tri_mem_cube;


//REALLY HACKY SOLUTION, SHOULD BE IN TRI_MEM_FUND, BUT THEN FROM_FILE AND TO_FILE BREAK :(
extern int * axis_sparse_to_normal;
extern size_t axis_sparse_len;


typedef struct tri_mem_fund {
  vert_index *    vert_to_index; //Array converts cube.points index to the index in our list
  arr3 *           vert_from_index; //Convert index to it's vertex
  vert_index **    sym_index;     //Cached symmetries
  int *            vert_fund_sym;  //Array that holds the symmetry number needed to convert this to fund domain
  size_t cube_len; //Total amount of points in the cube
  size_t fund_len; //Amount of points inside the fundamental domain
} tri_mem_fund;

typedef struct tri_mem_square {
  int dim_size; //Size of the data
} tri_mem_square;

typedef struct tri_mem_list
{
  unsigned char  *** t_arr; //Actual data array
  union {
    tri_mem_cube   mem_cube;
    tri_mem_tet    mem_tet;
    tri_mem_fund   mem_fund;
    tri_mem_square mem_square;
  };
  int dim;      //Original dimension
  int mode;     //Stores the type of mem_list (store every triangle in the cube, only in fundamental domain, or triangles in the unit tetrahedron)
  
} tri_mem_list;

#define MEM_LIST_CUBE 0
#define MEM_LIST_FUND 2
#define MEM_LIST_TET  3
#define MEM_LIST_CUBE_SPARSE 4
#define MEM_LIST_FUND_SPARSE 5
#define MEM_LIST_SQUARE 6

typedef vert_index tri_index[3];
typedef vert_index tet_index[4];

typedef struct tri_index_list
{
  tri_index * index_list;
  size_t len;
  int dim;
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
#define vertex_to_index_cube(vertex,dim) (vertex[0] * (dim+1)*(dim+1) + vertex[1] * (dim+1) + vertex[2])
#define vertex_to_index_cube_axis(vertex, axis, axis_len) (axis[vertex[0]] * axis_len * axis_len + axis[vertex[1]] * axis_len + axis[vertex[2]])
//Might want to make a seperate one for NON-sparse fund mem_list?
#define vertex_to_index_fund(vertex, vertex_index_array) (vertex_index_array[\
    vertex_to_index_cube_axis(vertex, axis_sparse_to_normal, axis_sparse_len)])

#define vertex_to_index_square(vert,p) (vert[0] * (p+1) + vert[1])
                                     
/* 
 * See indices_unique_*. This is a wrapper to convert vertices to a unique index directly
 */
                              
#define vertices_to_index_cube(v1, v2, v3, dim, indices) (\
          indices_unique_cube(vertex_to_index_cube(v1, dim),\
                              vertex_to_index_cube(v2, dim),\
                              vertex_to_index_cube(v3, dim),indices))

#define vertices_to_index_tet(v1, v2, v3, vert_to_index, indices) (\
          indices_unique_tet(vertex_to_index_tet(v1, vert_to_index),\
                             vertex_to_index_tet(v2, vert_to_index),\
                             vertex_to_index_tet(v3, vert_to_index),indices))

#define vertices_to_index_square(v1, v2, v3, p, indices) (\
          indices_unique_square(vertex_to_index_square(v1, p),\
                                vertex_to_index_square(v2, p),\
                                vertex_to_index_square(v3, p),indices))
                             
/*
 * These are wrappers to convert triangles (thus three vertices) to a unique index directly
 */
      
#define triangle_to_index_cube(triang, dim_mult, indices) (\
          vertices_to_index_cube(triang.vertices[0],triang.vertices[1], triang.vertices[2], dim_mult,indices))
          
#define triangle_to_index_tet(triang, vert_to_index, indices) (\
          vertices_to_index_tet( triang.vertices[0],triang.vertices[1], triang.vertices[2], vert_to_index,indices))

#define triangle_to_index_square(triang, p,indices) (\
    vertices_to_index_square(triang.vertices[0], triang.vertices[1], triang.vertices[2], p, indices))

/*
 * Returns whether this is a fundamental triangle (has point in fund domain, by using fund len)
 */
#define triangle_in_fund(idx1, idx2, idx3, fund_len) (((idx1 < fund_len) || (idx2 < fund_len) || (idx3 < fund_len)))

/*
 * Index - Vertex functions
 */
void indices_unique_square(vert_index idx1, vert_index idx2, vert_index idx3, tri_index indices);
void indices_unique_cube(vert_index idx1, vert_index idx2, vert_index idx3, tri_index indices);
void indices_unique_fund(vert_index idx1, vert_index idx2, vert_index idx3, tri_index indices);
void indices_unique_tet(vert_index idx1, vert_index idx2, vert_index idx3, tri_index indices);

/*
 * Functions that generate vertex -> index arrays
 */
void gen_axis(int dim, int sparse);
void gen_vertex_to_index_tet(int dim, vert_index_array * vertex_index, int *tet_len);


void vertex_from_index_cube(vert_index index, int dim,arr3 vertex) ;
triangle triangle_from_index_cube(tri_index indices, int dim);
triangle triangle_from_index_fund(tri_index indices,arr3 * index_vertex);

#define vertex_from_index_square(result, index, p) {\
  result[0] = (index) / (p+1);\
  result[1] = (index) % (p+1);}
triangle_2d triangle_from_index_square(tri_index indices, int p);


/*
 * mem_list init/free functions
 */
#define MEM_LIST_FALSE 0
#define MEM_LIST_TRUE 1

tri_mem_list mem_list_init(int dim, int mode, int init_value);
tri_mem_list mem_list_cube_init(int dim, int init_value, int sparse);
tri_mem_list mem_list_fund_init(int dim, int init_value, int mode);
tri_mem_list mem_list_tet_init(int dim, int init_value);
tri_mem_list mem_list_square_init(int p, int init_value);
void mem_list_free(tri_mem_list * list);



/*
 * mem_list operations
 */
int mem_list_dim_size(tri_mem_list * list, int dim, int idx1, int idx2);
int mem_list_row_empty(tri_mem_list * list,int i, int j);
size_t mem_list_count(tri_mem_list * list);
size_t mem_list_memory(tri_mem_list * list);
size_t mem_list_indices(tri_mem_list * list, tri_index_list * index_list);
triangle_list mem_list_to_triangle_list(tri_mem_list * list);

/*
 * mem_list conversion operators
 */
tri_mem_list mem_list_fund_to_cube(tri_mem_list * fund_list);

/*
 * mem_list file functions
 */
#define MEM_LIST_SAVE_FULL 0
#define MEM_LIST_SAVE_CLEAN 1
int mem_list_from_file(tri_mem_list * result, char * filename);
int mem_list_to_file(tri_mem_list * list, char * filename, int mode);









/*
 * mem_list set/get/clear +  cube functions
 */
void mem_list_cube_compress(tri_mem_list * list);
void mem_list_cube_clear(tri_mem_list * list, ptriangle triang);
void mem_list_cube_set(tri_mem_list * list, ptriangle triang);
int  mem_list_cube_get(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3);
int mem_list_cube_contains(tri_mem_list * list, ptriangle triang);
void mem_list_cube_clear_sym(tri_mem_list * list, ptriangle triang);
void mem_list_cube_set_sym(tri_mem_list * list, ptriangle triang);

/*
 * mem_list set/get/clear fund functions
 */
void mem_list_fund_clear(tri_mem_list * list, ptriangle triang);
void mem_list_fund_set(tri_mem_list * list, ptriangle triang);
int mem_list_fund_get(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3);
/*
 * mem_list set/get/clear tet functions
 */
int mem_list_tet_get(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3);
void mem_list_tet_set(tri_mem_list * list, ptriangle triang);
void mem_list_tet_clear(tri_mem_list * list, ptriangle triang);
/*
 * mem_list set/get/clear square functions
 */
int mem_list_square_get(tri_mem_list * list, arr2 v1, arr2 v2, arr2 v3);
void mem_list_square_set(tri_mem_list * list, ptriangle_2d triang);
void mem_list_square_clear(tri_mem_list * list, ptriangle_2d triang);
/*
 * old
 */
tri_mem_list mem_list_from_index_list(tri_index_list * list);
tri_mem_list mem_list_from_index_list_fund(tri_index_list * list);
int index_list_from_file(tri_index_list * result, char * filename);
int index_list_to_file(tri_index_list * list, char * filename);


#endif
