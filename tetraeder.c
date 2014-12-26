#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vector.h"
#include "triangle.h"
#include "tetraeder.h"
#include "omp.h"



int tet_tet_share_edge(ptetra t1, ptetra t2) {
  return (vert_vert_share_count(t1->vertices, 4, t2->vertices, 4) >= 2);
}

int tet_tet_share_facet(ptetra t1, ptetra t2){
  return (vert_vert_share_count(t1->vertices, 4, t2->vertices, 4) >= 3);
}

int tri_tet_share_edge(ptriangle tri, ptetra tet){
  return (vert_vert_share_count(tri->vertices,3, tet->vertices,4) >= 2);
}

int tri_tet_share_facet(ptriangle tri, ptetra tet) {
  return (vert_vert_share_count(tri->vertices,3, tet->vertices,4) == 3);
}

void tet_sides(ptetra tet, ptriangle triang) {
  arr3_to_triangle(tet->vertices[0], tet->vertices[1], tet->vertices[2], triang + 0);
  arr3_to_triangle(tet->vertices[0], tet->vertices[1], tet->vertices[3], triang + 1);
  arr3_to_triangle(tet->vertices[0], tet->vertices[2], tet->vertices[3], triang + 2);
  arr3_to_triangle(tet->vertices[1], tet->vertices[2], tet->vertices[3], triang + 3);
}
/*
 * Calculates the normal vector for each facet of tet. All vectors
 * point inwards or outwards.
 */
void tetra_normals(ptetra tet, arr3 * normals) {
  arr3 P[5]; //Edges
  subArr3(tet->vertices[1], tet->vertices[0], P[0]);
  subArr3(tet->vertices[2], tet->vertices[0], P[1]);
  subArr3(tet->vertices[3], tet->vertices[0], P[2]);
  
  subArr3(tet->vertices[2], tet->vertices[1], P[3]);
  subArr3(tet->vertices[3], tet->vertices[1], P[4]); 

  crossArr3(P[4],P[3], normals[0]); //Normal on facet 1,2,3 
  crossArr3(P[1],P[2], normals[1]); //Normal on facet 0,2,3
  crossArr3(P[2],P[0], normals[2]); //Normal on facet 0,1,3
  crossArr3(P[0],P[1], normals[3]); //Normal on facet 0,1,2
}
//The exact volume needs to be divided by six.

int tetra_volume(ptetra tet) {
  arr3 P[3]; //3 edge vectors from vertex a
  arr3 tmp;

  subArr3(tet->vertices[1], tet->vertices[0], P[0]); //b - a
  subArr3(tet->vertices[2], tet->vertices[0], P[1]); //c - a
  subArr3(tet->vertices[3], tet->vertices[0], P[2]); //d - a

  crossArr3(P[1],P[2], tmp);
  return abs(dotArr3(P[0], tmp));
}
/*
 * Function tests if the tetrahedron in the parameter is acute. Does this by calculating the normal
 * vectors and checking the dot products.
 */

int tetra_acute(ptetra tet) {
  arr3 P[5]; //Edges
  arr3 normals[4]; 
  subArr3(tet->vertices[1], tet->vertices[0], P[0]); //b - a
  subArr3(tet->vertices[2], tet->vertices[0], P[1]); //c - a
  subArr3(tet->vertices[3], tet->vertices[0], P[2]); //d - a
  subArr3(tet->vertices[2], tet->vertices[1], P[3]); //c - b
  subArr3(tet->vertices[3], tet->vertices[1], P[4]); //d - b
  
  crossArr3(P[4],P[3], normals[0]); //Normal on facet 1,2,3 
  crossArr3(P[1],P[2], normals[1]); //Normal on facet 0,2,3
  crossArr3(P[2],P[0], normals[2]); //Normal on facet 0,1,3
  crossArr3(P[0],P[1], normals[3]); //Normal on facet 0,1,2

  return (dotArr3(normals[0], normals[1]) < 0 &&
          dotArr3(normals[0], normals[2]) < 0 &&
          dotArr3(normals[0], normals[3]) < 0 &&
          dotArr3(normals[1], normals[2]) < 0 &&
          dotArr3(normals[1], normals[3]) < 0 &&
          dotArr3(normals[2], normals[3]) < 0);
}

/*
 * Does the same as the function above, but a optimized function.
 */
int tetra_acute_optimized(ptriangle tet, arr3 cube_pt) {
  arr3 P[5]; //Edges
  arr3 normals[4]; 
  /*
   * Calculate three edges of the tetrahedron
   */
  subArr3(tet->vertices[0], cube_pt, P[0]);
  subArr3(tet->vertices[1], cube_pt, P[1]);
  subArr3(tet->vertices[2], cube_pt, P[2]);
  
  /*
   * All Facets must have acute triangles, cheap test to rule this tetra out right now
   */
  if (dotArr3(P[0],P[1]) <= 0 || dotArr3(P[0], P[2]) <= 0 || dotArr3(P[1],P[2]) <= 0 )
    return 0;
  
  crossArr3(P[2],P[0], normals[2]); //Normal on facet 0,1,3
  crossArr3(P[0],P[1], normals[3]); //Normal on facet 0,1,2
  
  if (dotArr3(normals[2], normals[3]) >= 0)
    return 0;
    
  crossArr3(P[1],P[2], normals[1]); //Normal on facet 0,2,3
  if (dotArr3(normals[1], normals[2]) >= 0 ||
      dotArr3(normals[1], normals[3]) >= 0)
    return 0;
  
  subArr3(tet->vertices[1], tet->vertices[0], P[3]);
  subArr3(tet->vertices[2], tet->vertices[0], P[4]); 
  crossArr3(P[4],P[3], normals[0]); //Normal on facet 1,2,3 
  
  return (dotArr3(normals[0], normals[1]) < 0 &&
          dotArr3(normals[0], normals[2]) < 0 &&
          dotArr3(normals[0], normals[3]) < 0);
}


/*
 * Debug function
 */

void print_tetra(ptetra tet) {
  printf("[[%d,%d,%d],\n",tet->vertices[0][0],tet->vertices[0][1],tet->vertices[0][2]);
  printf(" [%d,%d,%d],\n",tet->vertices[1][0],tet->vertices[1][1],tet->vertices[1][2]);
  printf(" [%d,%d,%d],\n",tet->vertices[2][0],tet->vertices[2][1],tet->vertices[2][2]);
  printf(" [%d,%d,%d]]\n",tet->vertices[3][0],tet->vertices[3][1],tet->vertices[3][2]);
}

/*
 * Returns whether the facets of the tetrahedron, given by the triangle and apex, tetrahedron 
 * all ocur in the memory_list. 
 * 
 * Does not check wheter the facet given by triangle itself occurs in the memory_list. (Check before!).
 */

int facet_tetra_list(ptriangle triang, arr3 apex, data_list * data) {
  if (data->mode == DATA_MEM_LIST_FUND)
    return (mem_list_fund_get(&data->mem_list, triang->vertices[0], triang->vertices[1], apex) &&
            mem_list_fund_get(&data->mem_list, triang->vertices[0], triang->vertices[2], apex) &&
            mem_list_fund_get(&data->mem_list, triang->vertices[1], triang->vertices[2], apex));
  else if (data->mode == DATA_MEM_LIST_CUBE) 
    return (mem_list_cube_get(&data->mem_list, triang->vertices[0], triang->vertices[1], apex) &&
            mem_list_cube_get(&data->mem_list, triang->vertices[0], triang->vertices[2], apex) &&
            mem_list_cube_get(&data->mem_list, triang->vertices[1], triang->vertices[2], apex));
  else if (data->mode == DATA_MEM_LIST_TET)
    return (mem_list_tet_get(&data->mem_list, triang->vertices[0], triang->vertices[1],apex) &&
            mem_list_tet_get(&data->mem_list, triang->vertices[0], triang->vertices[2],apex) &&
            mem_list_tet_get(&data->mem_list, triang->vertices[1], triang->vertices[2],apex)); 
  else if (data->mode == DATA_TRI_LIST) 
    return (tri_list_get(&data->list, triang->vertices[0], triang->vertices[1],apex) &&
            tri_list_get(&data->list, triang->vertices[0], triang->vertices[2],apex) &&
            tri_list_get(&data->list, triang->vertices[1], triang->vertices[2],apex)); 
}

/*
 * Returns whether the facet given by triang is conform. i.e.
 * is a facet which is part of an acute tetrahedron of each "side" of it's plan, with a tetrahedron
 * having facets in the memory_list.
 * 
 */
int facet_conform(ptriangle triang, facet_acute_data * data) {
  if (data->store_acute_ind)
    data->acute_ind_len = 0;
  arr3 P[3];
  //Calculate the three edges of this triangle
  triangle_sides(triang->vertices[0], triang->vertices[1], triang->vertices[2],P);
  if (!triangle_P_acute(P))
    return 0;
  arr3 tri_normal;
  crossArr3(P[1], P[0], tri_normal); //Calculate normal on the triangle plane
  int  tri_d = dotArr3(tri_normal, triang->vertices[0]); //Find the constant specific for this plane
  
  //Calculate the vector perpendicular to each side and the normal. Normals on each side of the prism
  arr3 side_normals[3];
  arr3 side_d; //The constant expression for the plane equations of the sides
  //Convention, need to explain why it works?? Third point must be on the other side!!
  
  crossArr3(tri_normal, P[0], side_normals[0]);
  crossArr3(P[1], tri_normal, side_normals[1]);
  crossArr3(tri_normal, P[2], side_normals[2]);
  
  side_d[0] = dotArr3(side_normals[0], triang->vertices[0]); 
  side_d[1] = dotArr3(side_normals[1], triang->vertices[0]);
  side_d[2] = dotArr3(side_normals[2], triang->vertices[2]);
  
  
  data->boundary_triangle = data->boundary_func(triang, data->cube->dim); //Boundary plane only needs acute tetra on 1 side
  int dotprod;
  /*
   * Initalize to zero
   */
  data->acute_above = 0;
  data->acute_below = 0;
  /*
   * In case we work store the acute indices, we need to locally store
   * if we found above and below. As we later put them into the data struct.
   */
  int acute_above = 0;
  int acute_below = 0;
  for (unsigned short i = 0; i < data->cube->len; i++){
    dotprod = dotArr3(data->cube->points[i], tri_normal); //Calculate angle between normal and vector from plane to cube_pt
    //Check if current point lies on side of the triangle which isn't sharp yet and check if point lies in the
    //prism of possible sharp tetrahedron
    if (((dotprod > tri_d && !data->acute_above) ||  //Dotprod > tri_d implies that the point lies "above" the triangle (same side as normal) 
         (dotprod < tri_d && !data->acute_below)) && //Check if point lies in the prism
          dotArr3(data->cube->points[i],side_normals[0]) < side_d[0] &&
          dotArr3(data->cube->points[i],side_normals[1]) < side_d[1] &&
          dotArr3(data->cube->points[i],side_normals[2]) < side_d[2] &&    //Dotprod < tri_d implies lies below
          tetra_acute_optimized(triang,data->cube->points[i]) && //Tetrahedron is acute
          facet_tetra_list(triang, data->cube->points[i], data->data)) //All facets are in the conf_mem_list
      { //Passed all tests, we have found a correct tetrahedron on this side.

        if (data->store_acute_ind) { //We need to store all apices

	  data->acute_ind[data->acute_ind_len++] = i;
	  if (dotprod > tri_d)
	    acute_above = 1;
	  else
	    acute_below = 1;

        } else {

          if (dotprod > tri_d) 
            data->acute_above = 1;
          else 
            data->acute_below = 1;
          if ((data->acute_above && data->acute_below) || data->boundary_triangle)
            return 1;    
        }


      }
  }

  //If we stored all tetrahedrons, we now need to check if this facet is conform
  if (data->store_acute_ind) {
    data->acute_above = acute_above;
    data->acute_below = acute_below;

    if (data->acute_above && data->acute_below)
      return 1;

    if (data->boundary_triangle && (data->acute_above || data->acute_below))
      return 1;
  }

  return 0;  
}
//if save_file is set. The conf_mem_list is saved to this file every hour
#define save_interval 15*60

/*
 * Filters all non-conform facets in the conf_mem_list. Keeps looping
 * untill every facet is a conform facet, or the list is empty. This
 * turns the start set given by conf_mem_list into a conform set.
 * 
 * Saves the memory_list every save_interval, if save_file is set.
 */

void facets_conform_cube(data_list * data,int fund_domain, char * save_file){
  int changed = 1;
  tri_index indices;
  char tmp_file[100];
  if (save_file) 
    sprintf(tmp_file,"%s_tmp", save_file);
  
  
  size_t i,j,k;
  triangle cur_tri;
  tri_mem_list * conf_mem_list = &data->mem_list;
  facet_acute_data parameters;
  cube_points cube = gen_cube_points(conf_mem_list->dim);
  cube_points fund = gen_fund_points(conf_mem_list->dim);
  parameters.cube = &cube;
  parameters.store_acute_ind= 0;
  parameters.boundary_func = &triangle_boundary_cube;
  parameters.data = data;
  
  double time_start =0 , time_end = 0;
  while (changed) {
    changed = 0;
    printf("Starting conform loop with %zu facets.\n"
          , mem_list_count(conf_mem_list));   
    time_start = omp_get_wtime();

    if (!fund_domain) {
      #pragma omp parallel for  schedule(dynamic,conf_mem_list->dim) private(cur_tri, indices, i,j,k) firstprivate(parameters)
      for (i = 0; i < cube.len; i++) 
	for (j = 0; j < cube.len - i; j++)
	  if (conf_mem_list->t_arr[i][j])
	    for (k = 0; k < cube.len - j - i; k++)
	    {
	      indices[0] = i;
	      indices[1] = j;
	      indices[2] = k;
	      if (!GMI(conf_mem_list->t_arr, indices))
		continue; 

	      indices[1] = i + j;
	      indices[2] = i + j + k;
	      cur_tri = triangle_from_index_cube(indices, conf_mem_list->dim);

	      if (!facet_conform(&cur_tri,&parameters)) { //remove from conf_mem_list
		changed = 1;
		printf("What! Found a non-conform triangle :-(\n");
		mem_list_cube_clear(conf_mem_list, &cur_tri);
	      }
	  }
    } else {
      for (i = 0; i < fund.len; i++) {
	for (j = 0; j < cube.len; j++) {
	  for (k = j;k < cube.len; k++) //All combinations of three vertices in the tet
	  {
	    //Below is the same as indices_unique(i,j,k,indices) because i < j < k, already sorted

	    indices_unique_cube(vertex_to_index_cube(fund.points[i],conf_mem_list->dim),j,k,indices);
	    if (!GMI(conf_mem_list->t_arr,indices)) //Check if this index is still acute
	      continue;
	    cur_tri = (triangle) {{{fund.points[i][0],fund.points[i][1],fund.points[i][2]},
				  {cube.points[j][0],cube.points[j][1],cube.points[j][2]},
				  {cube.points[k][0],cube.points[k][1],cube.points[k][2]}}};
	    if (!facet_conform(&cur_tri,&parameters)) { //remove from list
	      changed = 1;
	      mem_list_cube_clear_sym(conf_mem_list, &cur_tri);
	    }
	  }
	}
      }
    }
    time_end   = omp_get_wtime();
    printf("Loop took %f seconds.\n\n", time_end-time_start);
  }
}

/*
 * See facets_conform_cube. Only this checks for all the triangles in the 
 * unit tetrahedron.
 */
void facets_conform_tet(data_list * data, char * save_file){
  int changed = 1;
  tri_index indices;
  char tmp_file[100];
  if (save_file) 
    sprintf(tmp_file,"%s_tmp", save_file);
  
  vert_index i,j,k;
  triangle cur_tri;
  tri_mem_list * conf_mem_list =&data->mem_list;
  facet_acute_data parameters;
  cube_points tet = gen_tet_points(conf_mem_list->dim);
  parameters.cube = &tet;
  parameters.boundary_func = &triangle_boundary_tet;
  parameters.data = data;
  parameters.store_acute_ind = 0;
  
  double time_start =0 , time_end = 0, time_save;
  while (changed) {
    changed = 0;
    printf("Starting conform loop with %zu facets.\n"
          , mem_list_count(conf_mem_list));   
    time_start = omp_get_wtime();
    time_save =  save_interval;
    #pragma omp parallel for schedule(dynamic) private(j,k,i,cur_tri,indices)  firstprivate(parameters)
    for (i = 0; i < tet.len; i++) {
      for (j = i + 1; j < tet.len; j++) {
        for (k = j + 1;k < tet.len; k++) //All combinations of three vertices in the tet
        {
          //printf("%d", omp_get_thread_num());
          //Below is the same as indices_unique(i,j,k,indices) because i < j < k, already sorted
          
          indices[0] = i;
          indices[1] = j - i - 1;
          indices[2] = k - j - 1;
          if (!GMI(conf_mem_list->t_arr,indices)) //Check if this index is still acute
            continue;
            
          cur_tri = (triangle) {{{tet.points[i][0],tet.points[i][1],tet.points[i][2]},
                                {tet.points[j][0],tet.points[j][1],tet.points[j][2]},
                                {tet.points[k][0],tet.points[k][1],tet.points[k][2]}}};
          if (!facet_conform(&cur_tri,&parameters)) { //remove from list
            changed = 1;
            CMI(conf_mem_list->t_arr,indices); 
          }
        }
      }
      if (save_file && //Do we want to save the file?
          omp_get_thread_num() == 0 && //Only let master save to the file
        ((omp_get_wtime() - time_start) > time_save)) { //Time to save current progress
        
        printf("Saving tmp file with +/- %zu facets.\n", mem_list_count(conf_mem_list));
        if (!mem_list_to_file(conf_mem_list, tmp_file, MEM_LIST_SAVE_CLEAN)) { //Save as tmp
          remove(tmp_file);
          time_save += save_interval / 2;
        } else {
          remove(save_file); //Remove the old progress file
          rename(tmp_file, save_file); //Rename tmp to new   
          time_save += save_interval;
        }
      }
    }
    time_end   = omp_get_wtime();
    printf("Conform loop took %f seconds.\n\n", time_end-time_start);
  }
  free(tet.points);
}

/*
 * See facets_conform_cube. Only stores triangles in the fundamental domain.
 */
void facets_conform_fund(data_list * data, char * save_file){
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
}

void facets_conform_tri_list(data_list * data, char * save_file){
  int changed = 1;
  char tmp_file[100];
  if (save_file) 
    sprintf(tmp_file,"%s_tmp", save_file);

  int l;
  size_t i,j,k;
  triangle cur_tri;
  tri_list * conf_list = &data->list;
  facet_acute_data parameters;
  cube_points cube = gen_cube_points(conf_list->dim);
  parameters.cube = &cube;
  parameters.boundary_func = &triangle_boundary_cube;
  parameters.data = data;
  parameters.store_acute_ind = 0;

  double time_start =0 , time_end = 0, time_save;
  while (changed) {
    changed = 0;
    printf("Starting conform loop with %zu facets.\n"
        , tri_list_count(conf_list));   
    time_start = omp_get_wtime();
    time_save = save_interval;

    #pragma omp parallel for schedule(dynamic,conf_list->dim) private(j,k,i,l,cur_tri)  firstprivate(parameters)
    for (i = 0; i < cube.len; i++) {
      for (j = i; j < cube.len; j++) {
        for (l = conf_list->t_arr[i][j-i].len - 1; l >= 0; l--) {  //Loop over all triangles (i,j,*)
          k = conf_list->t_arr[i][j-i].p_arr[l] +  j;
          //Just vertex_from_index_cube?
          cur_tri = (triangle) {{{cube.points[i][0],cube.points[i][1],cube.points[i][2]},
                                 {cube.points[j][0],cube.points[j][1],cube.points[j][2]},
                                 {cube.points[k][0],cube.points[k][1],cube.points[k][2]}}};
          if (!facet_conform(&cur_tri,&parameters)) {
            changed = 1;
            tri_list_remove(conf_list, &cur_tri, TRI_LIST_NO_RESIZE);
          }
        }
        if (omp_get_thread_num() == 0 && //Only let master save to the file
          ((omp_get_wtime() - time_start) > time_save))  //Time to save current progress
        {  
          printf("We have %zu facets left.\n", tri_list_count(conf_list));
          time_save += save_interval;
        }
      }
    }
    time_end   = omp_get_wtime();
    printf("Conform loop took %f seconds.\n\n", time_end-time_start);
  }
}
/*
 * Calls the correct conform function. Depens on the memory type
 */
void facets_conform(data_list * data, char * save_file){
  switch (data->mode) {
    case DATA_MEM_LIST_FUND:
      facets_conform_fund(data, save_file);
      break;
    case DATA_MEM_LIST_CUBE:
      facets_conform_cube(data,1, save_file);
      break;
    case DATA_MEM_LIST_TET:
      facets_conform_tet(data,save_file);
      break;
    case DATA_TRI_LIST:
      facets_conform_tri_list(data, save_file);
      break;
  }
}

/*
 * Counts the number of acute tetrahedra in [0,dim]^3.
 */
int tetrahedra_acute(int dim) {
  double time_start, time_end;
  time_start = omp_get_wtime();
  cube_points cube = gen_cube_points(dim);
  int c = 0;
  for (size_t j = 0; j < cube.len; j++)
    for (size_t k = j + 1; k < cube.len; k++)
      for (size_t l = k + 1; l < cube.len; l++)
        for (size_t m = l + 1; m < cube.len; m++) {
          tetra tet;
          copyArr3(tet.vertices[0], cube.points[j]);
          copyArr3(tet.vertices[1], cube.points[k]);
          copyArr3(tet.vertices[2], cube.points[l]);
          copyArr3(tet.vertices[3], cube.points[m]);
          if (tetra_acute(&tet))
            c++;
        }
  time_end = omp_get_wtime();
  printf("There are %d acute tetrahedra for p = %d\n", c, dim);
  printf("Calculating took %f seconds\n", time_end - time_start);
  return c;
}
