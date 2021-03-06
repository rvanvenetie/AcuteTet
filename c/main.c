#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <time.h>
#include "vector.h"
#include "triangle.h"
#include "tetraeder.h"
#include "mem_list.h"
#include "tri_list.h"
#include "omp.h"

/*
 * Loops from dimension 1 to 25 by default.
 * Parameters are "<type> <dim_1>" or "<type> <dim_1> <dim_2>".
 * Where type is either cube/fund/tet. Cube stores all the triangles of the cube, fund only
 * stores triangles in the fundamental domain and tet uses the unit tetrahedron instead of the cube.
 * If <dim_2> is also specified, it loops from dim_1 to dim_2.
 */
 
#define USE_FILE 1
#define REDIRECT_OUTPUT 0

#define CUBE_LOG         "output_cube_%d.log"
#define FUND_LOG         "output_fund_%d.log"
#define FUND_SPARSE_LOG  "output_fund_sparse_%d.log"
#define TET_LOG          "output_tet_%d.log"
#define TRI_LOG          "output_tri_%d.log"

#define FUND_SPARSE_TMP "/local/rvveneti/fund_sparse_conf_%d.tet"
#define FUND_TMP	"/local/rvveneti/fund_conf_%d.tet"
#define TET_TMP		"/local/rvveneti/tet_conf_%d.tet"
#define CUBE_TMP	"/local/rvveneti/cube_verz_%d.tet"
#define TRI_TMP		"/local/rvveneti/tri_conf_%d.tet"


#define FUND_SPARSE_DATA "data/fund_sparse_conf_%d.tet"
#define FUND_DATA	 "data/fund_conf_%d.tet"
#define TET_DATA	 "data/tet_conf_%d.tet"
#define CUBE_DATA	 "data/cube_conf_%d.tet"
#define TRI_DATA	 "data/tri_conf_%d.tet"


#define LOOP_FUND        DATA_MEM_LIST_FUND
#define LOOP_FUND_SPARSE DATA_MEM_LIST_FUND_SPARSE
#define LOOP_CUBE	 DATA_MEM_LIST_CUBE
#define LOOP_TET         DATA_MEM_LIST_TET
#define LOOP_TRI         DATA_TRI_LIST

int main(int argc, char *argv[]) {
  char tmp_file[70],log_file[70],data_file[70];
  data_list face_list;
  double time_start,time_end;
  
  
  int loop_type = LOOP_FUND; //By default loop over the fundamental domain
  int loop_start = 1; //Start with dimension 1
  int loop_end   = 25; //End with dimension 25
  int help = 0;
  if (argc == 3 || argc == 4) { //Check arguments
    if (argv[1][0] == 'c')
      loop_type = LOOP_CUBE;
    else if (argv[1][0] == 't')
      loop_type = LOOP_TET;
    else if (argv[1][0] == 's')
      loop_type = LOOP_FUND_SPARSE;
    else if (argv[1][0] == 'l')
      loop_type = LOOP_TRI;
    else if (!(argv[1][0] == 'f')) {
      help = 1;
    }
    loop_start = atoi(argv[2]);
    if (argc == 4)
      loop_end = atoi(argv[3]);
    else
      loop_end = loop_start;
  } else if (argc >= 2) 
    help = 1;
  if (help) {
      printf("Illegal parameters. \n");
      printf("Parameters are <type> <p> or <type> <start p> <end p>\n");
      printf("  <type> can be cube/fund/sparse/tet/list.\n");
      printf("   - Cube: Stores all the facets in the cube. \n");
      printf("   - Fund: Only stores facets with a point in the fundamental domain. \n");
      printf("   - Tet:  Tries to find a Conforme Verzameling for the unit tetrahedron instead of the cube. \n");
      printf("Example parameters \"fund 5 10\", this generates the Conforme Verzameling for p=5..10. \n");    
      exit(0);
  }
  printf("OpenMP information for this computer\n");
  printf("Dynamic threads:%d\n", omp_get_dynamic());
  printf("Maximum number of threads: %d\n",omp_get_max_threads());
  printf("Example parallel for loop: \n");
  //omp_set_num_threads(48);
  #pragma omp parallel for schedule(static,1)
  for (int i = 0; i < 70; i++)
    printf("(%d,%d)", omp_get_thread_num(),omp_get_num_threads());
  printf("\n");
   
  printf("Header tri_mem_list:%lu\n", sizeof(tri_mem_list));
  for (int dim = loop_start;dim < loop_end + 1; dim++) 
  {
      
    int i;
    if (loop_type == LOOP_FUND_SPARSE)
      i = 1 << dim;
    else
      i = dim;

    switch (loop_type) {
      case LOOP_FUND: 
        sprintf(log_file,  FUND_LOG,i); 
        sprintf(tmp_file,  FUND_TMP,i);         
        sprintf(data_file, FUND_DATA,i);         
        break;

      case LOOP_FUND_SPARSE:
        sprintf(log_file,  FUND_SPARSE_LOG,i); 
        sprintf(tmp_file,  FUND_SPARSE_TMP,i);         
        sprintf(data_file, FUND_SPARSE_DATA,i);         
        break;

      case LOOP_CUBE: 
        sprintf(log_file,  CUBE_LOG,i);
        sprintf(tmp_file,  CUBE_TMP,i);
        sprintf(data_file, CUBE_DATA,i);
        break;
      case LOOP_TET:
       sprintf(log_file,  TET_LOG,i);
       sprintf(tmp_file,  TET_TMP,i);
       sprintf(data_file, TET_DATA,i);
       break;
      case LOOP_TRI:
       sprintf(log_file,  TRI_LOG,i);
       sprintf(tmp_file,  TRI_TMP,i);
       sprintf(data_file, TRI_DATA,i);
       break;
    }
      
    if (REDIRECT_OUTPUT) {
      if (freopen(log_file,"a",stdout) == NULL)
        printf("Redirecting output failed\n");
      setvbuf(stdout, NULL,_IOLBF, 1024);
    }
    printf("\n\n\nGenerating the Conforme Verzameling for p=%d\n\n\n\n", i);
    
    time_start = omp_get_wtime();
    
    
    if (USE_FILE && data_list_from_file(&face_list, loop_type,tmp_file)){
      printf("Continuing previous data-set.\n");
    } else {
      printf("Initalizing new data-set.\n");
      face_list = data_list_init(i, loop_type, MEM_LIST_TRUE);
    }    
    data_list_to_file(&face_list, "data/tst123", MEM_LIST_SAVE_CLEAN);
    time_end   = omp_get_wtime();
    printf("Took %f seconds to init the memory list.\n", time_end - time_start);
    time_start = omp_get_wtime();
    printf("Size of the memory list is %zu bytes.\n", data_list_memory(&face_list));
    printf("Amount of start facets is %zu.\n", data_list_count(&face_list));
    printf("Start filtering triangles not conform:\n\n");
    printf("Cube_len %zu fund_len %zu\n", face_list.mem_list.mem_fund.cube_len, face_list.mem_list.mem_fund.fund_len);

    facets_conform(&face_list, tmp_file);

    time_end = omp_get_wtime();
    printf("\nAmount of conform facets after filtering: %zu\n", data_list_count(&face_list));    
    printf("Total calculation Conforme Verzameling took %f seconds\n\n",time_end - time_start);



    if (!data_list_to_file(&face_list, data_file, MEM_LIST_SAVE_CLEAN))
      printf("Failed to save conform data to file: %s\n", data_file);

    printf("Total memory used: %zu\n\n", data_list_memory(&face_list));

    data_list_free(&face_list);
  
    if (REDIRECT_OUTPUT)
      fclose(stdout);
  }

}
