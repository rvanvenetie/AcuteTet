#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <time.h>
#include "vector.h"
#include "triangle.h"
#include "tetraeder.h"
#include "mem_list.h"
#include "combinations.h"
#include "omp.h"
#include "triangulate.h"
/*
 * Loops from dimension 1 to 25 by default.
 * Parameters are "<type> <dim_1>" or "<type> <dim_1> <dim_2>".
 * Where type is either tet or fund (tetrahedron or fundamental region cube).
 * If <dim_2> is also specified, it loops from dim_1 to dim_2.
 */
 
#define USE_FILE 1
#define REDIRECT_OUTPUT 1

#define FUND_LOG "output_fund_%d.log"
#define TET_LOG "output_tet_%d.log"

#define FUND_TMP "/var/scratch/rvveneti/gez_verz_%d.tet"
#define TET_TMP  "/var/scratch/rvveneti/tet_verz_%d.tet"

#define FUND_DATA "data/gez_verz_%d.tet"
#define TET_DATA  "data/tet_verz_%d.tet"


int main(int argc, char *argv[]) {
  printf("Thread numbers:\n");
  #pragma omp parallel for
  for (int i = 0; i < 16; i++)
    printf("%d,", omp_get_thread_num());
  printf("\n");
  char filename[70];
  char log_file[70];
  /*
  ptriangulation opdeling = triangulate_cube_random(dim);
  while (!opdeling)
    opdeling = triangulate_cube_random(dim);
  exit(0);
  */
  tri_mem_list face_list;
  double time_start,time_end;
  
  int loop_tet = 1; //By default loop over the tetrahedron shape
  int loop_start = 1; //Start with dimension 1
  int loop_end   = 25; //End with dimension 25
  if (argc == 3 || argc == 4) { //Check arguments
    if (argv[1][0] == 'f')
      loop_tet = 0;
    else if (!(argv[1][0] == 't')) {
      printf("Illegal parameters. \n");
      exit(0);
    }
    loop_start = atoi(argv[2]);
    if (argc == 4)
      loop_end = atoi(argv[3]);
    else
      loop_end = loop_start;
  }
  
  for (int i = loop_start;i < loop_end + 1; i++) 
  {
    if (loop_tet)
      sprintf(log_file, TET_LOG, i);
    else
      sprintf(log_file, FUND_LOG,i);
      
    if (REDIRECT_OUTPUT) {
      if (freopen(log_file,"a",stdout) == NULL)
        printf("Redirecting output failed\n");
      setvbuf(stdout, NULL,_IOLBF, 1024);
    }
    printf("\n\n\nGenerating the gezellige verzameling for dimension %d\n\n\n\n", i);
    
    time_start = omp_get_wtime();
    
    if (loop_tet)
      sprintf(filename,TET_TMP,i);
    else
      sprintf(filename,FUND_TMP,i);
    
    
    if (USE_FILE && mem_list_from_file(&face_list, filename)) {
      printf("Continuing previous data-set.\n");
    } else {
      printf("Initalizing new data-set.\n");
      if (loop_tet)
        face_list = mem_list_init_tet(i, MEM_LIST_TRUE);
      else 
        face_list = mem_list_init_fund2(i, MEM_LIST_TRUE);
    }    
    time_end   = omp_get_wtime();
    printf("Took %f seconds to init the memory list.\n", time_end - time_start);
    time_start = omp_get_wtime();
    printf("Size of the memory list is %zu bytes.\n", mem_list_memory(&face_list));
    printf("Start filtering triangles not acute or not gezellig.\n\n");
    

    facets_face2face(&face_list, filename);
      
    time_end = omp_get_wtime();
    printf("\nAmount of sharp facets after face2face filter: %zu\n", mem_list_count(&face_list));    
    printf("Total calculation gezellige verzameling took %f seconds\n\n",time_end - time_start);
    
    if (loop_tet)
      sprintf(filename,TET_DATA,i);
    else
      sprintf(filename,FUND_DATA,i);
    
    
    if (!mem_list_to_file(&face_list, filename, MEM_LIST_SAVE_CLEAN))
      printf("Failed to save face2face data to file: %s\n", filename);
    //printf("Total memory used: %zu\n\n", mem_list_memory(&face_list));
    
    mem_list_free(&face_list);
    if (REDIRECT_OUTPUT)
      fclose(stdout);
  }
}
