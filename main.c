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


#define USE_FILE 1
#define FACE2FACE 0
#define REDIRECT_OUTPUT 1
#define LOOP



void test_multithread(void) {
  //omp_set_num_threads(4);
  printf( "Starting execution with %d threads:\n",
  omp_get_num_threads ());
  #pragma omp parallel
  {
    #pragma omp single
    printf("Thread %d\n", omp_get_thread_num()); 
    #pragma omp for
    for(int i = 0; i < 10; i++)
    {
      printf("Thread %d has number %d\n",omp_get_thread_num(), i);
    }
  }
}

int main(int argc, char *argv[]) {
  //test_multithread();
  //omp_set_num_threads(4);
  //exit(00);
  char filename[70];
  char log_file[70];
  arr3 dim =  {35,30,30};
  size_t len;
  /*
  ptriangulation opdeling = triangulate_cube_random(dim);
  while (!opdeling)
    opdeling = triangulate_cube_random(dim);
  exit(0);
  */
  //unsigned short *** vertex_index;
  //cube_points tet = gen_tet_points(dim[0], &vertex_index);
  //printf("%zu\n", tet.len);
 // exit(0);
  tri_mem_list fund_list, face_list;
  //triangle_list tri_list;
  tri_index_list idx_list;
  //tetra_list tet_list;
  
  double time_start,time_end;
 // tet_list = acute_tetrahedra_recur(dim);
 // exit(0);
  int i;
  int c = 0;
  #ifdef LOOP
  for (i = 23;i < 25; i++, c++) 
  #endif
  {
    if (argc == 2) {
      if (c== 1)
        exit(0); 
      i = atoi(argv[1]);
    }
    sprintf(log_file, "output_tet_%d.log",i);
    if (REDIRECT_OUTPUT) {
      freopen(log_file,"a",stdout);
      setvbuf(stdout, NULL,_IOLBF, 1024);
    }
    printf("\n\n\nGenerating the gezellige verzameling for dimension %d\n\n\n\n", i);
    
    time_start = omp_get_wtime();
    sprintf(filename,"/var/scratch/rvveneti/tet_verz_%d.tet",i);
    if (USE_FILE && mem_list_from_file(&face_list, filename)) {
      printf("Continuing previous data-set.\n");
    } else {
      printf("Initalizing new data-set.\n");
      face_list = mem_list_init_tet(i, MEM_LIST_TRUE); 
    }    
    time_end   = omp_get_wtime();
    printf("Took %f seconds to init the memory list.\n", time_end - time_start);
    time_start = omp_get_wtime();
    printf("Size of the memory list is %zu bytes.\n", mem_list_count(&face_list));
    printf("Start filtering triangles not acute or not gezellig.\n\n");
    facets_face2face_tet(&face_list, filename); //Start filtering
    time_end = omp_get_wtime();
    printf("\nAmount of sharp facets after face2face filter: %zu\n", mem_list_count(&face_list));    
    printf("Total calculation gezellige verzameling took %f seconds\n\n",time_end - time_start);
    
    sprintf(filename,"data/tet_verz_%d.tet",i);
    printf("Memory before clean up: %zu\n", mem_list_memory(&face_list));
    if (!mem_list_to_file(&face_list, filename, MEM_LIST_SAVE_CLEAN))
      printf("Failed to save face2face data to file: %s\n", filename);
    printf("Total memory used: %zu\n\n", mem_list_memory(&face_list));
    
    mem_list_free(&face_list);
    fclose(stdout);
    continue;
    
    
    printf("Dimension: %d. Gathering fundamental data.\n", i);
    sprintf(filename,"/var/scratch/rvveneti/fund_data_%d.tet",i);
    if (USE_FILE && mem_list_from_file(&fund_list, filename))
      printf("Succesfully data from file%s\n", filename);
    else {
      time_start = omp_get_wtime();
      fund_list = facets_cube_acute(i);
      time_end   = omp_get_wtime();
      printf("Calculated data in %f seconds\n", time_end - time_start);  
      if (!mem_list_to_file(&fund_list, filename,MEM_LIST_SAVE_CLEAN)) 
        printf("Failed to save to file %s\n", filename);      
    }
    printf("Amount of sharp facets in fundamental domain:%zu\n", mem_list_count(&fund_list));
    printf("Memory used by fundamental data: %zu\n", mem_list_memory(&fund_list));
    if (!FACE2FACE) {
      #ifdef LOOP
      mem_list_free(&fund_list);
      continue;
      #else
      exit(0);
      #endif
    }
    sprintf(filename,"data/face_data_%d.tet",i);
    printf("\nFiltering fund data with the face2face property\n");
    if (USE_FILE && mem_list_from_file(&face_list,filename )) {
      mem_list_free(&fund_list);
      printf("Loaded face2face data from file: %s\n", filename);
    } else{ 
      face_list = fund_list;
      time_start = omp_get_wtime();
      facets_face2face_fund(&face_list,NULL);
      time_end   = omp_get_wtime();
      printf("Calculated data in %f seconds\n", time_end   - time_start);  
      printf("Memory before clean up: %zu\n", mem_list_memory(&face_list));
      mem_list_clean(&face_list);
      if (!mem_list_to_file(&face_list, filename, MEM_LIST_SAVE_CLEAN))
        printf("Failed to save face2face data to file: %s\n", filename);
    }
    printf("Amount of sharp facets after face2face filter: %zu\n", mem_list_count(&face_list));    
    printf("Total memory used: %zu\n\n", mem_list_memory(&face_list));
    mem_list_free(&face_list);
  }
}
