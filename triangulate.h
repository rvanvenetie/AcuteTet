#ifndef TRIANGULATE_H
#define TRIANGULATE_H
#include "limits.h"
#include "triangle.h"
#include "tetraeder.h"

typedef struct triangulation {
  ptriangle boundaries;
  ptetra    tetraeders;
  size_t    len_bound, len_tetra;  
  arr3 dim;
} triangulation, ptriangulation;

int tetra_tetra_disjoint(ptetra t1, ptetra t2);

#endif
