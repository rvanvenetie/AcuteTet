#ifndef EDGE_H
#define EDGE_H

#include "vector.h"
#include "triangle.h"
#include "edge_list.h"

typedef struct triangulation {
	pedge      bound_edge;
	int 		   bound_len;
	
	ptriangle_2d  tri;
	int        tri_len;
	int p;
} triangulation, * ptriangulation;

typedef struct {
	edge_matrix * mat;
  int acute_above, acute_below;
	int boundary_edge;

	/*
	 * Below is needed if you want to store the actual acute triangles
	 */
	int store_acute_ind;
  vert_index * acute_ind;
  int acute_ind_len;
} edge_conform_parameters;
#endif
