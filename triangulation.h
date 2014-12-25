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

#define TOUCH 2
#define DISJOINT 1
#define INTERSECT 0

int tet_tet_disjoint(ptetra t1, ptetra t2);
int tri_tet_disjoint(ptriangle tri, ptetra tet);
void triangulation_free(ptriangulation triang);
void triangulation_print(ptriangulation triang);
int triangulation_from_file(ptriangulation triang, char * filename);
int triangulation_to_file(ptriangulation triang, char * filename);

#define TRIANGULATE_NO_CONFORM 0
#define TRIANGULATE_CONFORM 1
triangulation triangulate_cube(data_list * data,  char * tmp_triang_file, char * tmp_data_file);
#endif
