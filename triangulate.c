#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#include "vector.h"
#include "triangle.h"
#include "tetraeder.h"
#include "edge_list.h"
#include "mem_list.h"
#include "tri_list.h"
#include "triangulation.h"
#include "omp.h"

#define REDIRECT_OUTPUT 1
#define LOG "triang_%d_%s.log"
#define TRIANG_FILE "/local/rvveneti/triangulation"
#define TRIANG_TMP_FILE "/local/rvveneti/triang_tmp"
#define TRIANG_TET_TMP_FILE "/local/rvveneti/triangulation_triangles_tmp.tet"
#define TRIANG_TET_FILE "/local/rvveneti/triangulation_triangles.tet"
int main(int argc, char *argv[]) {
  int p = 31;
  tri_list t_list;

  tri_list_from_file(&t_list,"/local/rvveneti/tri_list_31.tet");
  edge_matrix mat = edge_matrix_init(p, 0);
  int row_size = edge_matrix_row_size(&mat);
  /*
   * All triangles in x = 0 plane are given by the indices
   * i < j < k < row_size. Note that an edge (i,j) exists if len > 0 for tri[i][j-i].len
   */
  //Loop over all edges in x = 0 plane
  size_t totalcnt = 0;
  size_t edge_cnt = 0;
	for (int i = 0; i < row_size; i++) 
		for (int j = 0; j < row_size - i; j++) {
      for (int k = 0; k < t_list.t_arr[i][j].len; k++) {
        tri_index idx;
        idx[0] = i;
        idx[1] = i + j;
        idx[2] = i + j + t_list.t_arr[i][j].p_arr[k];
        if (idx[2] < row_size) { //Third point also in the x = 0 plane
          totalcnt++;
          //Edge_count, increase if mat.val = 0
          edge_cnt += 1 - mat.val[idx[0]][idx[1]] + 1 - mat.val[idx[1]][idx[2]] + 1 - mat.val[idx[0]][idx[2]];

          //Set all edges to one
          mat.val[idx[0]][idx[1]] = mat.val[idx[1]][idx[0]] = 1;
          mat.val[idx[1]][idx[2]] = mat.val[idx[2]][idx[1]] = 1;
          mat.val[idx[0]][idx[2]] = mat.val[idx[2]][idx[0]] = 1;
        }
      }
    }
  triangle cur_tri;
	//Loop over all triangles in the square
  size_t tri_cnt = 0;
  size_t tri_cnt_bla = 0;
	for (int i = 0; i < row_size; i++) {
    vertex_from_index_cube(i,mat.p, cur_tri.vertices[0]);
		for (int j = i; j < row_size; j++) {
      vertex_from_index_cube(j,mat.p, cur_tri.vertices[1]);
			for (int k = j; k < row_size; k++) {
        vertex_from_index_cube(k,mat.p, cur_tri.vertices[2]);
				//If triangle is acute and the sides are contained in the matrix
				if (arr3_triangle_acute(cur_tri.vertices[0],cur_tri.vertices[1],cur_tri.vertices[2]) && 
						mat.val[i][j] &&
						mat.val[i][k] &&
						mat.val[j][k])
				{
          if (tri_list_contains(&t_list, &cur_tri))
            tri_cnt_bla++;
          tri_cnt++;
				}

			}
		}
	}
  printf("Tricnt    %zu\n", tri_cnt_bla);
  printf("Tri cnt   %zu\n", tri_cnt);
  printf("Totalcnt= %zu\n", totalcnt);
  printf("Edgecnt = %zu\n", edge_cnt);
  printf("Edge_mat= %zu\n", edge_matrix_count(&mat));
  edge_matrix_cosy_count(&mat, &totalcnt, &edge_cnt);
  printf("Cosy tri =%zu\n", totalcnt);
  printf("Cosy ara =%zu\n", edge_cnt);
  edge_matrix_to_file(&mat, "/var/scratch/rvveneti/plane.mat");
  return 1;
	/*
		 tri_mem_list fund_list;
		 tri_list     t_list;
		 fprintf(stderr,"Start loading\n");
		 if (mem_list_from_file(&fund_list, "/var/scratch/rvveneti/fund_conf_31.tet"))
		 fprintf(stderr,"Loaded succesfuly\n");
		 else
		 {
		 fprintf(stderr,"Failed loading..\n");
		 return 0;
		 }
		 tri_mem_list cube_list;
		 cube_list = mem_list_fund_to_cube(&fund_list);
		 if (mem_list_to_file(&cube_list, "/local/rvveneti/cube_conf_31.tet", MEM_LIST_SAVE_CLEAN))
		 fprintf(stderr,"Stored succesfully\n");
		 else
		 fprintf(stderr, "Failed storing\n");
		 return 0;
		 t_list = mem_list_to_tri_list(&fund_list);
		 fprintf(stderr,"Not failed in mem_list_to_tri_list!\n");
		 if (tri_list_to_file(&t_list, "/local/rvveneti/tri_list_31.tet"))
		 fprintf(stderr,"Stored succesfully\n");
		 else
		 fprintf(stderr,"Failed storing..  tri list\n");
		 return 0;
	 */
	srand(1234);
  srand(time(NULL));
	char log_file[100];
	int dim = 10;
	data_list list;
	triangulation triang;
	while (1) {
		if (argc == 2 && isdigit(argv[1][0]))
			dim = atoi(argv[1]);
    int mode = DATA_TRI_LIST;
		if (argc == 2 && !isdigit(argv[1][0]))
		{
      if (strstr(argv[1], "cube")) //Cube in the name, must be cube
        mode = DATA_MEM_LIST_CUBE;
			fprintf(stderr,"Loading from: %s\n", argv[1]);
			if (data_list_from_file(&list, mode, argv[1]))
				fprintf(stderr,"Loaded succesfully from file\n");
			else {
				fprintf(stderr,"Failed loading from file\n");
				return -1;
			}
			dim = data_list_dim(&list);
		} else {
			fprintf(stderr,"Creating new tri_list\n");
			list = data_list_init(dim, mode, MEM_LIST_TRUE);
		}
    if (mode == DATA_TRI_LIST)
      sprintf(log_file, LOG, dim, "tri_list");
    else 
      sprintf(log_file, LOG, dim, "cube");
		if (REDIRECT_OUTPUT) {
			if (freopen(log_file,"a",stdout) == NULL)
				printf("Redirecting output failed\n");
			setvbuf(stdout, NULL,_IOLBF, 1024);
		}
		printf("Triangulation for p = %d\n", dim);
		printf("Amount of triangles in the list = %zu\n", data_list_count(&list));
		printf("Memory of the triangle list = %zu\n",     data_list_memory(&list));

		printf("Finding triangulation.\n");
		//triang = triangulate_cube_random(&list);
    do {
      triang = triangulate_cube(&list,  TRIANG_TMP_FILE, TRIANG_TET_TMP_FILE);
      if (triang.bound_len == 0) //No boundary triangles!
      {
        printf("Triangulation found!");
        triangulation_to_file(&triang, TRIANG_FILE);
        data_list_to_file(&list, TRIANG_TET_FILE, MEM_LIST_SAVE_CLEAN);
        break;
      } 

      triangulation_free(&triang);
    } while (mode == DATA_TRI_LIST && tri_list_update_from_file(&list.list, argv[1]));
		data_list_free(&list);
    fprintf(stderr, "Restarting loop\n\n");
	}
	triangulation_free(&triang);
	data_list_free(&list);
	return 0;
}
