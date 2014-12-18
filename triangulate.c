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

#define REDIRECT_OUTPUT 1
#define LOG "triang_%d.log"
#define TRIANG_FILE "/local/rvveneti/triangulation"
#define TRIANG_TMP_FILE "/local/rvveneti/triang_tmp"
#define TRIANG_TET_TMP_FILE "/local/rvveneti/triangulation_triangles_tmp.tet"
#define TRIANG_TET_FILE "/local/rvveneti/triangulation_triangles.tet"
int main(int argc, char *argv[]) {
  char log_file[100];
  int dim = 10;
  tri_mem_list list;
  ptriangulation triang;
  while (1) {
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
    //printf("Amount of triangles in the list = %zu\n", mem_list_count(&list));
    printf("Memory of the triangle list = %zu\n", mem_list_memory(&list));

    triang = triangulate_cube(&list, TRIANG_TMP_FILE, TRIANG_TET_TMP_FILE);
    if (triang->bound_len == 0) //No boundary triangles!
    {
      printf("Triangulation found!");
      triangulation_to_file(triang, TRIANG_FILE);
      mem_list_to_file(&list, TRIANG_TET_FILE, MEM_LIST_SAVE_CLEAN);
      break;
    } 

    triangulation_free(triang);
    mem_list_free(&list);
  }
  triangulation_free(triang);
  mem_list_free(&list);
  return 0;
}
