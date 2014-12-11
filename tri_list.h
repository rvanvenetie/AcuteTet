#ifndef TRI_LIST_H 
#define TRI_LIST_H 

#include "mem_list.h"

typedef struct int_arr{
  unsigned short * p_arr;
  int len;
} int_arr;
/*
 * This is a storage class for triangles. It does so by creating a 3D-array.
 * The first two axis indicate the first two point of the triangles. The last 
 * array is the array of all the triangles having these first points.

 * The last array is sorted, this way we can quickly look up if a triangle
 * is contained in this list.
 */
typedef struct tri_list
{
  int_arr **  t_arr; //The actual data
  int dim;      //Original dimension
} tri_list;

#define TRI_LIST 1337

int int_arr_index(int_arr * arr, int i);
tri_list tri_list_init(int dim, int init_value);
void tri_list_insert(tri_list * list, ptriangle  triang);
int tri_list_contains(tri_list * list, ptriangle  triang);
int tri_list_get(tri_list * list, arr3 v1, arr3 v2, arr3 v3);
size_t tri_list_count(tri_list * list);
void tri_list_remove(tri_list * list, ptriangle  triang);
void tri_list_free(tri_list * list);

int tri_list_from_file(tri_list * result, char * filename);
int tri_list_to_file(tri_list * list, char * filename);
tri_list mem_list_to_tri_list(tri_mem_list * list);
#endif
