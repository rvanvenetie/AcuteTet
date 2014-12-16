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

  /*
   * In case we want to reconstruct all the conform tetrahedrons, we can store the indices,
   * in cube points, of the apices that form the acute tetrahedron with this base triangle.
   * One must ensure that the array allows enough storage. (upper bound given by the amount of points
   * in cube).
   */
  int store_acute_ind; 
  vert_index * acute_ind;
  int acute_ind_len;
} facet_acute_data;


int tet_tet_share_edge(ptetra t1, ptetra t2);
int tet_tet_share_facet(ptetra t1, ptetra t2);
int tri_tet_share_edge(ptriangle tri, ptetra tet);
int tri_tet_share_facet(ptriangle tri, ptetra tet);

void tet_sides(ptetra tet, ptriangle triang);
void tetra_normals(ptetra tet, arr3 * normals);
int tetra_acute(ptetra tet);
int tetra_acute_optimized(ptriangle tet, arr3 cube_pt);
int tetrahedra_acute(int dim);
void facets_conform(data_list * conform_list, char * save_file);
int facet_conform(ptriangle triang, facet_acute_data * data);
void print_tetra(ptetra tet);
#endif
