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

#define REDIRECT_OUTPUT 0
#define LOG "triang_%d.log"
int main(int argc, char *argv[]) {
  char log_file[100];
  int dim = 10;
  tri_mem_list list;
  if (argc == 2 && isdigit(argv[1][0]))
    dim = atoi(argv[1]);
  if (argc == 2 && !isdigit(argv[1][0]))
  {
    printf("Loading from: %s\n", argv[1]);
    if (mem_list_from_file(&list,argv[1]))
      printf("Loaded succesfully from file\n");
    else {
      printf("Failed loading from file\n");
      return -1;
    }
  } else {
    printf("Creating new tri_list\n");
    list = mem_list_init(dim, MEM_LIST_CUBE_SPARSE, MEM_LIST_TRUE);
  }
  sprintf(log_file, LOG, list.dim);
  if (REDIRECT_OUTPUT) {
    if (freopen(log_file,"a",stdout) == NULL)
      printf("Redirecting output failed\n");
    setvbuf(stdout, NULL,_IOLBF, 1024);
  }
  printf("Triangulation for p = %d\n", list.dim);
  printf("Amount of triangles in the list = %zu\n", mem_list_count(&list));
  printf("Memory of the triangle list = %zu\n", mem_list_memory(&list));

  ptriangulation triang;
  triang = triangulate_cube(&list);
  if (triang) {
    printf("Triangulation found!");
    triangulation_free(triang);
  }

  mem_list_free(&list);
  return 0;
}
