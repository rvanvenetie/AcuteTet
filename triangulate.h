#ifndef TRIANGULATE_H
#define TRIANGULATE_H
#include "limits.h"
#include "triangle.h"
#include "tetraeder.h"

typedef struct triangulation {
  ptriangle bound_tri;
  size_t    bound_len;
  ptetra    tetra;
  size_t    tetra_len;  

  
  int dim;
} triangulation, * ptriangulation;

#define DISJOINT 1
#define INTERSECT 0

void triangulation_free(ptriangulation triang);
int tetra_tetra_disjoint(ptetra t1, ptetra t2);
ptriangulation triangulate_cube_random(int dim);
#endif
