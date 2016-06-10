#ifndef EDGE_LIST_H
#define EDGE_LIST_H
#include "vector.h"
#include "mem_list.h"

typedef vert_index edge_index[2]; //[0] = x, [1] = y

typedef struct edge{
	arr2 vertices[2];
} *pedge;

typedef struct {
	unsigned char ** val; //2D Data array
	int p;
} edge_matrix;

#ifdef INLINE_MACROS
  #define xy_to_index_square(x,y,p) (x * (p+1) + y)
  #define edge_to_index_square(result, edge, p) {\
	  result[0] = vertex_to_index_square(edge->vertices[0], p);\
	  result[1] = vertex_to_index_square(edge->vertices[1], p);\
	}
	#define edge_from_index_square(result, indices, p) {\
		vertex_from_index_square(result->vertices[0], indices[0], p);\
		vertex_from_index_square(result->vertices[1], indices[1], p);\
	}
#else
vert_index vertex_to_index_square(arr2 vert, int p);
void edge_to_index_square(edge_index result, pedge edge, int p);
void  vertex_from_index_square(arr2 result, vert_index index, int p);
void edge_from_index_square(pedge result, edge_index indices, int p);
#endif
int edge_matrix_row_size(edge_matrix * mat);


edge_matrix edge_matrix_init(int p, unsigned char init_val);
edge_matrix edge_matrix_copy(edge_matrix * mat);
void edge_matrix_free(edge_matrix * mat);


void edge_matrix_set(edge_matrix * mat, pedge edge);
void edge_matrix_clear(edge_matrix * mat, pedge edge);
void edge_matrix_clear_sym(edge_matrix * mat, pedge edge);
int  edge_matrix_contains(edge_matrix * mat, pedge edge);
unsigned char edge_matrix_get(edge_matrix * mat, arr2 v1, arr2 v2);
size_t edge_matrix_count(edge_matrix * mat );

void edge_matrix_cosy_count(edge_matrix * mat, size_t * tri_cnt, size_t * total_area_x2);

int edge_matrix_from_file(edge_matrix * result, char  * filename) ;
int edge_matrix_to_file(edge_matrix * mat, char * filename) ;
#endif
