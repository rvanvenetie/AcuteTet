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

int tet_tet_disjoint(ptetra t1, ptetra t2);
int tri_tet_disjoint(ptriangle tri, ptetra tet);
void triangulation_free(ptriangulation triang);
#endif
