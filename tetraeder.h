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
  tri_mem_list * acute_list;
  
 //Store the actual tetrahedrons acute above and acute below
  ptetra tetra_above;
  int    tetra_above_len;
  ptetra tetra_below;
  int    tetra_below_len;
  
  //Store whether facet has one tetra acute above, acute below
  int acute_above, acute_below;
  //Store whether this triangle lies on a boundary (thus only is acute_above or acute_below)
  int boundary_triangle;
} facet_acute_data;


void tetra_normals(ptetra tet, arr3 * normals);
int tetra_acute(ptetra tetra);
tetra_list acute_tetrahedra(arr3 dim);
tetra_list acute_tetrahedra_recur(arr3 dim);
tri_mem_list acute_triangles_tetra(arr3 dim);
triangle_list acute_triangles_tetra_old(arr3 dim);
int triangle_tetra_acute(ptriangle triang, cube_points * cube, ptriang_tetra_result res  , tri_mem_list *acute_list);
void print_tetra(ptetra tet);
void mem_list_face2face_old(tri_mem_list * acute_list, int sym);
void mem_list_face2face(tri_mem_list * acute_list);


//Different modes
#define FACET_ACUTE 0
#define FACET_ACUTE_LIST 1
#define FACET_ACUTE_TETRA 2


void facets_face2face(tri_mem_list * acute_list, char * save_file);
tri_mem_list facets_cube_acute(int dim);
int facet_cube_acute(ptriangle triang, facet_acute_data * data, int mode);
#endif
