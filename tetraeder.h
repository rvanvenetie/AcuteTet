#ifndef TETRAEDER_H
#define TETRAEDER_H

#include "vector.h"
#include "triangle.h"
#include "stdlib.h"
#include "mem_list.h"
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

/*
 * Struct holds whether this triangle has atleast one acute tetrahedron "above"
 * and "below". Where above/below is defined by the normal on this triangle.
 * Furthermore holds whether this triangle is on the boundary of the cube.
 */
typedef struct
{
  int acute_above, acute_below;
  int boundary_above, boundary_below;
} triang_tetra_result, *ptriang_tetra_result;



typedef struct facet_acute_data {
  cube_points * cube;
  tri_mem_list * conform_list;
 
  //Store whether facet has one tetra acute above, acute below
  int acute_above, acute_below;
  
  //Store whether this triangle lies on a boundary (thus only is acute_above or acute_below)
  int (*boundary_func)(ptriangle, int);
  //int (*facets_list)(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3);

  int boundary_triangle;
} facet_acute_data;


void tetra_normals(ptetra tet, arr3 * normals);
int tetra_acute(ptetra tet);
int tetrahedra_acute(int dim);
void facets_conform(tri_mem_list * conform_list, char * save_file);
#endif
