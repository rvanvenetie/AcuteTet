#ifndef TRI_LIST_H 
#define TRI_LIST_H 

#include "vector.h"
#include "mem_list.h"


typedef struct int_arr_old{
  unsigned short * p_arr;
  size_t len;
} int_arr_old;

typedef struct tri_list_old
{
  int_arr_old **  t_arr; //The actual data
  int dim;      //Original dimension
} tri_list_old;
/*
 * This is a storage class for triangles. It does so by creating a 3D-array.
 * The first two axis indicate the first two point of the triangles. The last 
 * array is the array of all the triangles having these first points.

 * The last array is sorted, this way we can quickly look up if a triangle
 * is contained in this list.
 */

typedef struct int_arr{
  vert_index * p_arr; //Array of points
  unsigned short len;     //Amount of points we currently hold
  unsigned short data_len;//Amount of points we can store
} int_arr, * p_int_arr;

typedef struct tri_list
{
  int_arr ** t_arr;
  int dim;
} tri_list;
#define TRI_LIST 1337

int int_arr_index(int_arr * arr, int i);
tri_list tri_list_init(int dim, int init_value);
int tri_list_insert(tri_list * list, ptriangle  triang, int resize);
int tri_list_contains(tri_list * list, ptriangle  triang);
int tri_list_get(tri_list * list, arr3 v1, arr3 v2, arr3 v3);
size_t tri_list_count(tri_list * list);
size_t tri_list_memory(tri_list * list);
int tri_list_dim_size(tri_list* list, int axis, int idx1, int idx2);
#define TRI_LIST_NO_RESIZE 0
#define TRI_LIST_RESIZE 1
void tri_list_remove(tri_list * list, ptriangle  triang, int resize);
void tri_list_free(tri_list * list);
void tri_list_empty(tri_list * list);
void tri_list_resize(tri_list * list);

void tri_list_validate(tri_list * list);

int tri_list_from_file(tri_list * result, char * filename);
int tri_list_to_file(tri_list * list, char * filename);
tri_list mem_list_to_tri_list(tri_mem_list * list);

/*
 * This should be moved to a different unit..
 * Want to make one `mother' union that allows for all different data types
 * to be called by the same function.
 */

#define DATA_MEM_LIST_FUND 0
#define DATA_MEM_LIST_TET  1
#define DATA_MEM_LIST_CUBE 2
#define DATA_TRI_LIST      3

typedef struct data_list{
  int mode;
  union {
    tri_mem_list  mem_list;
    tri_list  list; 
  };
} data_list;

data_list data_list_init(int dim, int mode, int init_value);
int data_list_from_file(data_list * result, int mode, char * filename);
int data_list_to_file(data_list * list, char * filename, int mode);
size_t data_list_count(data_list * list);
size_t data_list_memory(data_list * list);
int data_list_dim(data_list * list);
int data_list_dim_size(data_list * list, int dim, int idx1, int idx2);

void data_list_free(data_list * list);
#endif

