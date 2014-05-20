#ifndef TRIANGULATE_H
#define TRIANGULATE_H
#include "limits.h"
#include "triangle.h"
#include "tetraeder.h"

typedef struct triangulation {
  ptriangle boundaries;
  ptetra    tetraeders;
  size_t    bound_len, tetra_len;  
  arr3 dim;
  arr2 dim_mult;
} triangulation, * ptriangulation;


void triangulation_free(ptriangulation triang);
int tetra_tetra_disjoint(ptetra t1, ptetra t2);
ptriangulation triangulate_cube_random(arr3 dim);
#endif
