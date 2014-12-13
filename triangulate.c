#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#include "vector.h"
#include "triangle.h"
#include "tetraeder.h"
#include "mem_list.h"
#include "tri_list.h"
#include "triangulation.h"
#include "omp.h"

int main(int argc, char *argv[]) {
  int dim = 10;
  tri_list list;
  if (argc == 2 && isdigit(argv[1][0]))
    dim = atoi(argv[1]);
  if (argc == 2 && !isdigit(argv[1][0]))
  {
    printf("Loading from: %s\n", argv[1]);
    if (tri_list_from_file(&list,argv[1]))
      printf("Loaded succesfully from file\n");
    else {
      printf("Failed loading from file\n");
      return -1;
    }
  } else {
    printf("Creating new tri_list\n");
    list = tri_list_init(dim, MEM_LIST_TRUE);
  }
  printf("Triangulation for p = %d\n", list.dim);
  printf("Amount of triangles in the list = %zu\n", tri_list_count(&list));
  printf("Memory of the triangle list = %zu\n", tri_list_memory(&list));

  ptriangulation triang;
  triang = triangulate_cube(&list);
  if (triang) {
    printf("Triangulation found!");
    triangulation_free(triang);
  }


  return 0;
}
