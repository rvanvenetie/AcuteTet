#include "conform.h"

/*
 * This class generates a conform set for a given cube
 */

void ConformCube::fund() {
  /*
  int changed = 1;
  tri_index indices;
  char tmp_file[100];
  if (save_file) 
    sprintf(tmp_file,"%s_tmp", save_file);


  tri_mem_list * conf_mem_list = &data->mem_list;
  arr3 * vert_from_index = conf_mem_list->mem_fund.vert_from_index;         
  size_t i,j,k;
  triangle cur_tri;
  facet_acute_data parameters;
  cube_points cube;
  if (conf_mem_list->mode == MEM_LIST_FUND)
    cube = gen_cube_points(conf_mem_list->dim);
  else if (conf_mem_list->mode == MEM_LIST_FUND_SPARSE)
    cube = gen_cube_sparse_points(conf_mem_list->dim);

  size_t fund_len = conf_mem_list->mem_fund.fund_len; //Cache
  parameters.cube = &cube;
  parameters.boundary_func = &triangle_boundary_cube;
  parameters.data = data;
  parameters.store_acute_ind = 0;

  double time_start =0 , time_end = 0, time_save;
  while (changed) {
    changed = 0;
    printf("Starting conform loop with %zu facets.\n"
        , mem_list_count(conf_mem_list));   
    time_start = omp_get_wtime();
    time_save = save_interval;

#pragma omp parallel for schedule(dynamic) private(j,k,i,cur_tri,indices)  firstprivate(parameters)
    for (i = 0; i < fund_len; i++) {
      for (j = i + 1; j < cube.len; j++) {
        for (k = j + 1; k < cube.len; k++) //All combinations of three vertices in the tet
        {
          //Below is the same as indices_unique(i,j,k,indices) because i < j < k, already sorted

          indices[0] = i;
          indices[1] = j - i - 1;
          indices[2] = k - j - 1;
          if (!GMI(conf_mem_list->t_arr,indices)) //Check if this index is still acute
            continue;

          cur_tri = (triangle) {{{vert_from_index[i][0],vert_from_index[i][1],vert_from_index[i][2]},
            {vert_from_index[j][0],vert_from_index[j][1],vert_from_index[j][2]},
            {vert_from_index[k][0],vert_from_index[k][1],vert_from_index[k][2]}}};

          if (!facet_conform(&cur_tri,&parameters)) { //remove from list
            changed = 1;
            mem_list_fund_clear(conf_mem_list, &cur_tri);
          }
        }
        if (save_file && //Do we want to save the file?
            omp_get_thread_num() == 0 && //Only let master save to the file
            ((omp_get_wtime() - time_start) > time_save))  //Time to save current progress
        {  
          printf("Saving tmp file with +/- %zu facets.\n", mem_list_count(conf_mem_list));
          remove(save_file); //Remove the old progress file
          if (!mem_list_to_file(conf_mem_list, save_file, MEM_LIST_SAVE_CLEAN)) { //Save as tmp
            remove(save_file);
            time_save += save_interval / 2;
          } else {
            time_save += save_interval;
          }
        }
      }
    }
    time_end   = omp_get_wtime();
    printf("Conform loop took %f seconds.\n\n", time_end-time_start);
  }
  free(cube.points);

  */
}
