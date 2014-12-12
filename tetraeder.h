#ifndef TETRAEDER_H
#define TETRAEDER_H

#include "vector.h"
#include "triangle.h"
#include "stdlib.h"
#include "mem_list.h"
#include "tri_list.h"
/*
 * Tetraeder has 4 vertices, each consists of 3 coordinates.
 */
typedef struct tetra
{
    arr3 vertices[4];
} tetra,*ptetra;

/*
 * Struct holds a list of all acute tetraeders for a given dimension
 */
typedef struct 
{
  ptetra t_arr;
  size_t len;
  arr3 dim;   
} tetra_list;



typedef struct facet_acute_data {
  data_list *  data;
  cube_points * cube;
  //Store whether facet has one tetra acute above, acute below
  int acute_above, acute_below;
  
  //Store whether this triangle lies on a boundary (thus only is acute_above or acute_below)
  int (*boundary_func)(ptriangle, int);
  //int (*facets_list)(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3);

  int boundary_triangle;
  
  int store_tetra;
  ptetra tet_above;
  ptetra tet_below;
  size_t tet_above_len;
  size_t tet_below_len;
} facet_acute_data;


void tetra_normals(ptetra tet, arr3 * normals);
int tetra_acute(ptetra tet);
int tetra_acute_optimized(ptriangle tet, arr3 cube_pt);
int tetrahedra_acute(int dim);
void facets_conform(data_list * conform_list, char * save_file);
int facet_conform(ptriangle triang, facet_acute_data * data);
void print_tetra(ptetra tet);
#endif
