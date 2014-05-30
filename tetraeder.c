#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "combinations.h"
#include "vector.h"
#include "triangle.h"
#include "tetraeder.h"
#include "omp.h"

/*
 * All normals point inward or outward... PROOF?
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
  
  /*
  // * Enforce outward pointing normals
  arr3 test_vec;
  subArr3(tet->vertices[0], tet->vertices[1],test_vec);
  if (dotArr3(normals[0],test_vec) > 0) {//Inward point
    negArr3(normals[0]);
    negArr3(normals[1]);
    negArr3(normals[2]);
    negArr3(normals[3]);
  }
  */
}

/*
 * Function tests if the tetrahedron in the parameter is acute. Calculates
 * the six normals on the faces. Then it calculates the dot product between
 * each pair and checks if the angle is acute.
 */

int tetra_acute(ptetra tet) {
  arr3 P[5]; //Edges
  arr3 normals[4]; 
  /*
   * Calculate three edges of the tetrahedron
   */
  subArr3(tet->vertices[1], tet->vertices[0], P[0]);
  subArr3(tet->vertices[2], tet->vertices[0], P[1]);
  subArr3(tet->vertices[3], tet->vertices[0], P[2]);
  
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
  
  subArr3(tet->vertices[2], tet->vertices[1], P[3]);
  subArr3(tet->vertices[3], tet->vertices[1], P[4]); 
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
 * Returns whether the facets of this tetrahedron all ocur in the memory_list (ie they must
 * be acute and have an acute tetra on each side).
 * 
 * The tetrahedron is made out of the base facet triang and the apex new_vertex. Function
 * assumes that the triang is in the memory_list (check before!)
 */

int facet_tetra_list(ptriangle triang, arr3 new_vertex, tri_mem_list * acute_list) {
  if (acute_list->mode == MEM_LIST_FUND2)
    return (mem_list_get_fund2(acute_list, triang->vertices[0], triang->vertices[1], new_vertex) &&
            mem_list_get_fund2(acute_list, triang->vertices[0], triang->vertices[2], new_vertex) &&
            mem_list_get_fund2(acute_list, triang->vertices[1], triang->vertices[2], new_vertex));
  else if (acute_list->mode == MEM_LIST_TET)
    return (mem_list_get_tet(acute_list, triang->vertices[0], triang->vertices[1],new_vertex) &&
            mem_list_get_tet(acute_list, triang->vertices[0], triang->vertices[2],new_vertex) &&
            mem_list_get_tet(acute_list, triang->vertices[1], triang->vertices[2],new_vertex));    
  else
    return (mem_list_get_fund(acute_list, triang->vertices[0], triang->vertices[1], new_vertex) &&
            mem_list_get_fund(acute_list, triang->vertices[0], triang->vertices[2], new_vertex) &&
            mem_list_get_fund(acute_list, triang->vertices[1], triang->vertices[2], new_vertex));    

}

/*
 * Add tetra to a tetra array
 */
void tetra_add_array(tetra tetra_to_add, ptetra  * tetra_array, int * len) {
  (*len)++;
  *tetra_array = (ptetra) realloc(*tetra_array, (*len) * sizeof(tetra));
  if (*tetra_array == NULL) {
    puts("Error allocating memory!!");
    exit(1);
  }
  (*tetra_array)[*len - 1] = tetra_to_add;      
}

/*
 * Returns whether the facet given by triang is acute. A acute facet is defined
 * to be a facet which is part of an acute tetrahedron on each "side" of it's plane. (Boundary
 * planes only need an acute tetrahedron on the inside). 
 * 
 * Different modes:
 *   -FACET_ACUTE:   Just checks if this facet has at least one acute tetrahedron on each side.
 *   -FACET_ACUTE_TETRA: Returns two arrays with the acute tetrahedron on each side
 *      of the facet.
 *   -FACET_ACUTE_LIST:  Only allows acute tetrahedron for which each of it's
 *     facets are part of the memory_list
 * 
 */
int facet_cube_acute(ptriangle triang, facet_acute_data * data, int mode) {
  /*
   * Every facet of an acute tetrahedron needs to be acute. If this facet is not even acute
   * we may directly stop checking this facet as it's never going to be part of an acute tetrahedron
   */
  if (!mat3_triangle_acute(triang->vertices)) 
    return 0;
    
  mat3 P;
  //Calculate the three edges of this triangle
  triangle_sides(triang->vertices[0], triang->vertices[1], triang->vertices[2],P);
  arr3 tri_normal;
  crossArr3(P[1], P[0], tri_normal); //Calculate normal on the triangle plane
  int  tri_d = dotArr3(tri_normal, triang->vertices[0]); //Find the constant specific for this plane
  
  //Calculate the vector perpendicular to each side and the normal. Normals on each side of the prism
  mat3 side_normals;
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
  data->acute_above = 0;
  data->acute_below = 0;
  data->tetra_above_len = 0;
  data->tetra_above = NULL;
  data->tetra_below_len = 0;
  data->tetra_below = NULL;
  
  for (size_t i = 0; i < data->cube->len; i++){
    dotprod = dotArr3(data->cube->points[i], tri_normal); //Calculate angle between normal and vector from plane to cube_pt
    //Check if current point lies on side of the triangle which isn't sharp yet and check if point lies in the
    //prism of possible sharp tetrahedron
    
    if (((dotprod > tri_d && !data->acute_above) ||  //Dotprod > tri_d implies that the point lies "above" the triangle (same side as normal) 
         (dotprod < tri_d && !data->acute_below))     //Dotprod < tri_d implies lies below
          && //Check if point lies in the prism
          dotArr3(data->cube->points[i],side_normals[0]) < side_d[0] &&
          dotArr3(data->cube->points[i],side_normals[1]) < side_d[1] &&
          dotArr3(data->cube->points[i],side_normals[2]) < side_d[2]) 
      {
      tetra test_tetra;
      memcpy(test_tetra.vertices, &data->cube->points[i], sizeof(arr3));
      memcpy(test_tetra.vertices + 1, triang->vertices, 3 * sizeof(arr3)); 
      if (tetra_acute(&test_tetra)) {
        //All the facets must be in the acute_list
        if ((mode == FACET_ACUTE_LIST) && !facet_tetra_list(triang, data->cube->points[i], data->acute_list))
          continue;   
        
        //Explicitly create a list of the acute tetrahedron
        if (mode == FACET_ACUTE_TETRA) {
          if (dotprod > tri_d)
            tetra_add_array(test_tetra, &data->tetra_above, &data->tetra_above_len);
          else
            tetra_add_array(test_tetra, &data->tetra_below, &data->tetra_below_len);
        }
        //We only need to know if tetrahedron above and below acute
        else {
          if (dotprod > tri_d) 
            data->acute_above = 1;
          else 
            data->acute_below = 1;
          if ((data->acute_above && data->acute_below) || data->boundary_triangle)
            return 1;
        }
      }
    }
  }
  if (mode == FACET_ACUTE_TETRA) {
    data->acute_above = (data->tetra_above_len > 0);
    data->acute_below = (data->tetra_below_len > 0);
    if ((data->acute_above && data->acute_below) || (data->boundary_triangle && (data->acute_above || data->acute_below)))
      return 1;
  }
  return 0;  
}

/*
 * Returns a memory list with all the acute facets (see above for definition) in the cube
 * with grid 0..dim
 */
tri_mem_list facets_cube_acute(int dim) {
  tri_mem_list result = mem_list_init_fund(dim,MEM_LIST_FALSE);
  arr3 tmpdim = {dim,dim,dim};
  cube_points fund = gen_fund_points(dim);
  cube_points cube = gen_cube_points(tmpdim);
  size_t comb_len;
  Dindex * comb = combinations_list(cube.len, 2, &comb_len);
  
  size_t j,k,i,l;
  triangle cur_tri;
  facet_acute_data parameters;
  parameters.cube = &cube;
  parameters.boundary_func = &triangle_boundary_cube;
  #pragma omp parallel for schedule(guided) private(l,j,k,i,cur_tri) firstprivate(parameters)

  for (l = 0; l < comb_len; l++) {//Loop through all combinations
    j = comb[l*2];
    k = comb[l*2+1];
    for (i = 0; i < fund.len; i++){    
      cur_tri = (triangle) {{{fund.points[i][0],fund.points[i][1],fund.points[i][2]},
                             {cube.points[j][0],cube.points[j][1],cube.points[j][2]},
                             {cube.points[k][0],cube.points[k][1],cube.points[k][2]}}};
                      
      if (facet_cube_acute(&cur_tri,&parameters,FACET_ACUTE))
        mem_list_set_fund(&result, &cur_tri);
    }
  }
  free(comb);
  free(cube.points);
  free(fund.points);
  return result;
}

//if save_file is set. The acute_list is saved to this file every hour
#define save_interval 60*60

/*
 * Filters the memory_list for face-to-face property. It loops over every set
 * triangle in the fundamental domain of the cube and checks:
 * 1) Whether there is an acute tetrahedron on each side of this facet
 * 2) Whether this tetrahedron has facets that are also in this list
 * 
 * This should result in a list of facets which have an acute tetrahedron
 * on each side of the facet that is also build out of facets having this property.
 * i.e. face-to-face.
 * 
 * Saves the memory_list every save_interval, if save_file is set.
 */
void facets_face2face_fund(tri_mem_list * acute_list, char * save_file){
  cube_points fund = gen_fund_points(acute_list->dim[0]);
  cube_points cube = gen_cube_points(acute_list->dim);
  int changed = 1;
  size_t i;
  tri_index indices;
  char tmp_file[100];
  if (save_file) 
    sprintf(tmp_file,"%s_tmp", save_file);
          
  vert_index j, k;
  triangle cur_tri;
  facet_acute_data parameters;
  parameters.cube = &cube;
  parameters.boundary_func = &triangle_boundary_cube;
  parameters.acute_list = acute_list;
  
  double time_start =0 , time_end = 0, time_save = save_interval;
  while (changed) {
    changed = 0;
    printf("Face2face loop with %zu acute triangles from thread %d.\n"
          , mem_list_count(acute_list), omp_get_thread_num());
    time_start = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic) private(j,k,i,cur_tri, indices)  firstprivate(parameters)
    for (j = 0; j < cube.len; j++) {
      for (k = j+1; k < cube.len; k++) {
        for (i = 0; i < fund.len; i++){ //Against all fundamental points
          indices_unique_fund(vertex_to_index_cube(fund.points[i],acute_list->mem_fund.dim_mult), j,k, acute_list, indices); 
          if (!GMI(acute_list->t_arr,indices)) //Check if this index is still acute
            continue;
          cur_tri = (triangle) {{{fund.points[i][0],fund.points[i][1],fund.points[i][2]},
                               {cube.points[j][0],cube.points[j][1],cube.points[j][2]},
                               {cube.points[k][0],cube.points[k][1],cube.points[k][2]}}};
          if (!facet_cube_acute(&cur_tri,&parameters,FACET_ACUTE_LIST)) { //remove from list
            changed = 1;
            mem_list_clear_fund(acute_list, &cur_tri);
          }
        }
        if (save_file && //Do we want to save the file?
            omp_get_thread_num() == 0 && //Only let master save to the file
          ((omp_get_wtime() - time_start) > time_save)) { //Time to save current progress
          
          printf("Saving tmp file with +/- %zu facets.\n", mem_list_count(acute_list));
          if (!mem_list_to_file(acute_list, tmp_file, MEM_LIST_SAVE_CLEAN)) { //Save as tmp
            remove(tmp_file);
            time_save += save_interval / 2;
          } else {
            remove(save_file); //Remove the old progress file
            rename(tmp_file, save_file); //Rename tmp to new   
            time_save += save_interval;
          }
        }
      } 
    }
    time_end   = omp_get_wtime();
    printf("Loop took %f seconds.\n\n", time_end-time_start);
  }
  
  free(cube.points);
  free(fund.points);
}

/*
 * See facets_face2face_fund. Only this checks for all the triangles in the 
 * unit tetrahedron.
 */
void facets_face2face_tet(tri_mem_list * acute_list, char * save_file){
  int changed = 1;
  tri_index indices;
  char tmp_file[100];
  if (save_file) 
    sprintf(tmp_file,"%s_tmp", save_file);
  
  vert_index i,j,k;
  triangle cur_tri;
  facet_acute_data parameters;
  cube_points tet = gen_tet_points(acute_list->dim[0]);
  parameters.cube = &tet;
  parameters.boundary_func = &triangle_boundary_tet;
  parameters.acute_list = acute_list;
  
  
  double time_start =0 , time_end = 0, time_save = save_interval;
  while (changed) {
    changed = 0;
    printf("Face2face loop with %zu acute triangles from thread %d.\n"
          , mem_list_count(acute_list), omp_get_thread_num());
    time_start = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic) private(j,k,i,cur_tri,indices)  firstprivate(parameters)
    for (i = 0; i < tet.len; i++) {
      for (j = i + 1; j < tet.len; j++) {
        for (k = j + 1;k < tet.len; k++) //All combinations of three vertices in the tet
        {
          //printf("%d", omp_get_thread_num());
          //Below is the same as indices_unique_tet(i,j,k,indices) because i < j < k, already sorted
          
          indices[0] = i;
          indices[1] = j - i - 1;
          indices[2] = k - j - 1;
          if (!GMI(acute_list->t_arr,indices)) //Check if this index is still acute
            continue;
            
          cur_tri = (triangle) {{{tet.points[i][0],tet.points[i][1],tet.points[i][2]},
                                {tet.points[j][0],tet.points[j][1],tet.points[j][2]},
                                {tet.points[k][0],tet.points[k][1],tet.points[k][2]}}};
          if (!facet_cube_acute(&cur_tri,&parameters,FACET_ACUTE_LIST)) { //remove from list
            changed = 1;
            CMI(acute_list->t_arr,indices); 
          }
        }
      }
      if (save_file && //Do we want to save the file?
          omp_get_thread_num() == 0 && //Only let master save to the file
        ((omp_get_wtime() - time_start) > time_save)) { //Time to save current progress
        
        printf("Saving tmp file with +/- %zu facets.\n", mem_list_count(acute_list));
        if (!mem_list_to_file(acute_list, tmp_file, MEM_LIST_SAVE_CLEAN)) { //Save as tmp
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
    printf("Loop took %f seconds.\n\n", time_end-time_start);
  }
  free(tet.points);
}

/*
 * See face2face tet/fund
 */
void facets_face2face_fund2(tri_mem_list * acute_list, char * save_file){
  int changed = 1;
  tri_index indices;
  char tmp_file[100];
  if (save_file) 
    sprintf(tmp_file,"%s_tmp", save_file);
  
  
  arr3 * vert_from_index = acute_list->mem_fund2.vert_from_index;         
  size_t i,j,k;
  triangle cur_tri;
  facet_acute_data parameters;
  cube_points cube = gen_cube_points(acute_list->dim);
  size_t fund_len = acute_list->mem_fund2.fund_len; //Cache
  parameters.cube = &cube;
  parameters.boundary_func = &triangle_boundary_cube;
  parameters.acute_list = acute_list;
  
  double time_start =0 , time_end = 0, time_save = save_interval;
  while (changed) {
    changed = 0;
    printf("Face2face loop with %zu acute triangles from thread %d.\n"
          , mem_list_count(acute_list), omp_get_thread_num());
    time_start = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic) private(j,k,i,cur_tri,indices)  firstprivate(parameters)
    for (i = 0; i < fund_len; i++) {
      for (j = i + 1; j < cube.len; j++) {
        for (k = j + 1;k < cube.len; k++) //All combinations of three vertices in the tet
        {
          //Below is the same as indices_unique_tet(i,j,k,indices) because i < j < k, already sorted
          
          indices[0] = i;
          indices[1] = j - i - 1;
          indices[2] = k - j - 1;
          if (!GMI(acute_list->t_arr,indices)) //Check if this index is still acute
            continue;
            
          cur_tri = (triangle) {{{vert_from_index[i][0],vert_from_index[i][1],vert_from_index[i][2]},
                                 {vert_from_index[j][0],vert_from_index[j][1],vert_from_index[j][2]},
                                 {vert_from_index[k][0],vert_from_index[k][1],vert_from_index[k][2]}}};
          if (!facet_cube_acute(&cur_tri,&parameters,FACET_ACUTE_LIST)) { //remove from list
            changed = 1;
            mem_list_clear_fund2(acute_list, &cur_tri);
          }
        }
        if (save_file && //Do we want to save the file?
            omp_get_thread_num() == 0 && //Only let master save to the file
          ((omp_get_wtime() - time_start) > time_save))  //Time to save current progress
        {  
          printf("Saving tmp file with +/- %zu facets.\n", mem_list_count(acute_list));
          if (!mem_list_to_file(acute_list, tmp_file, MEM_LIST_SAVE_CLEAN)) { //Save as tmp
            remove(tmp_file);
            time_save += save_interval / 2;
          } else {
            remove(save_file); //Remove the old progress file
            rename(tmp_file, save_file); //Rename tmp to new   
            time_save += save_interval;
          }
        }
      }
    }
    time_end   = omp_get_wtime();
    printf("Loop took %f seconds.\n\n", time_end-time_start);
  }
}

void facets_face2face(tri_mem_list * acute_list, char * save_file){
  if (acute_list->mode == MEM_LIST_FUND)
    facets_face2face_fund(acute_list, save_file);
  else if (acute_list->mode == MEM_LIST_FUND2)
    facets_face2face_fund2(acute_list, save_file);
  else if (acute_list->mode == MEM_LIST_TET)
    facets_face2face_tet(acute_list,save_file);  
}
/*
 * Generate a list of acute tetrahedron for dim given as argument
 */
tetra_list acute_tetrahedra(arr3 dim){
  tetra_list result = {NULL, 0, {dim[0],dim[1],dim[2]}};
  size_t count = 0;
  ptetra t_arr = NULL;
  size_t len = 0;
  arr3 *cube_pts;
  cube_points cube = gen_cube_points(dim);
  len = cube.len;
  cube_pts = cube.points;
  
  printf("Amount of cube_pts: %zu. Total amount of combinations checked: %zu\n", len, len*len*len*len);
  Dindex	*arr;
  arr = revdoor_init(len, 4);  
  if (arr == 0) {
    puts("Error involving revdoor_init");
    exit(1);
    return result;
  }
  //Loop through all combinations
  do {
    size_t i = arr[0], j = arr[1], k = arr[2], l = arr[3];  
    tetra cur_tet = (tetra) {{{cube_pts[i][0],cube_pts[i][1],cube_pts[i][2]},
                              {cube_pts[j][0],cube_pts[j][1],cube_pts[j][2]},
                              {cube_pts[k][0],cube_pts[k][1],cube_pts[k][2]},
                              {cube_pts[l][0],cube_pts[l][1],cube_pts[l][2]}}};
    if (tetra_acute(&cur_tet)) {
      count++;
      t_arr = (ptetra) realloc(t_arr, count * sizeof(tetra));
      if (t_arr == NULL) {
        puts("Error allocating memory!!");
        exit(1);
      }
      t_arr[count - 1] = cur_tet;
    }
  } while (revdoor_next(arr));
  revdoor_free(arr);
  
    
  result.t_arr = t_arr;
  result.len = count;
  free(cube_pts);
  return result;
}

/*
 * Does the same as above, but works by recursion (less acute checks)
 */
tetra_list acute_tetrahedra_recur(arr3 dim){
  tetra_list result = {NULL, 0, {dim[0],dim[1],dim[2]}};
  tetra_list prev;
  int maxdim, axis;
  size_t count = 0;
  ptetra t_arr = NULL;
  size_t len = 0;
  arr3 *cube_pts;
  maxdim = maxArr3(dim, &axis);
  if (maxdim == 0 || (dim[0] + 1) * (dim[1] + 1) * (dim[2]+1) < 4)
    return result;
  arr3 new_dim = {dim[0],dim[1],dim[2]};
  new_dim[axis]--;
  prev = acute_tetrahedra_recur(new_dim);
   
  cube_points cube = gen_cube_points(dim);
  len = cube.len;
  cube_pts = cube.points;
  
	Dindex	*arr;
  arr = revdoor_init(len, 4);  
  if (arr == 0) {
    puts("Error involving revdoor_init");
    exit(1);
    return result;
  }
  //Loop through all combinations
  do {
    size_t i = arr[0], j = arr[1], k = arr[2], l = arr[3];
    //Find all acute tetrahedra with two vertices on both "edge-planes"
    
    /*
     * It is possible to skip this if:
     *   --first vertex has point with axis = 0
     *   --second vertex has point with axis = maxdim
     *   --last two vertices must be combination of all points, except points on
     *       axis = 0 with lower value than first vertex
     *       and except axis = maxdim with lower value than second vertex
     * 
     *   Need to define "lower". This is to exclude double combinations. 
     * 
     *  Possible to skip cube_pts.
     *  x coordinate is index div (dim[1]*dim[2])
     *  y coordinate is 
     *  z coordinate is index mod dim[2]??
     */
    if ((cube_pts[i][axis] == 0 || cube_pts[j][axis] == 0 || 
       cube_pts[k][axis] == 0 || cube_pts[l][axis] == 0) &&
      (cube_pts[i][axis] == maxdim || cube_pts[j][axis] == maxdim ||
       cube_pts[k][axis] == maxdim || cube_pts[l][axis] == maxdim))
    {
      tetra cur_tet = (tetra) {{{cube_pts[i][0],cube_pts[i][1],cube_pts[i][2]},
                                {cube_pts[j][0],cube_pts[j][1],cube_pts[j][2]},
                                {cube_pts[k][0],cube_pts[k][1],cube_pts[k][2]},
                                {cube_pts[l][0],cube_pts[l][1],cube_pts[l][2]}}};
      if (tetra_acute(&cur_tet)) {
        count++;
        t_arr = (ptetra) realloc(t_arr, count * sizeof(tetra));
        if (t_arr == NULL) {
          puts("Error allocating memory!!");
          exit(1);
        }
        t_arr[count - 1] = cur_tet;
      }
    }
        
  } while (revdoor_next(arr));
  revdoor_free(arr);  
  
  
  //Gather all old tetrahedra with point on expanding edge-plane
  ptetra prev_edges = NULL;
  size_t edge_count = 0;
  for (size_t i = 0; i < prev.len; i++) {
    if (prev.t_arr[i].vertices[0][axis] == maxdim - 1 ||
        prev.t_arr[i].vertices[1][axis] == maxdim - 1 ||
        prev.t_arr[i].vertices[2][axis] == maxdim - 1 ||
        prev.t_arr[i].vertices[3][axis] == maxdim - 1) //Tetrahedra was on expanding edge-plane
    {
      edge_count++;
      prev_edges = (ptetra) realloc(prev_edges, edge_count * sizeof(tetra));
      if (prev_edges == NULL) {
        puts("Error allocating memory!!");
        exit(1);
      }
      prev_edges[edge_count - 1] = prev.t_arr[i];
      //Shift this edge tetrahedra
      prev_edges[edge_count - 1].vertices[0][axis] += 1;
      prev_edges[edge_count - 1].vertices[1][axis] += 1;
      prev_edges[edge_count - 1].vertices[2][axis] += 1;
      prev_edges[edge_count - 1].vertices[3][axis] += 1;
    }
  }
  //Now the total set of acute tetrahedra constis of the old set + old_edge shifted + new edge calculated above  
  t_arr = realloc(t_arr,(count + edge_count + prev.len) * sizeof(tetra));
  
  memcpy(t_arr + count , prev_edges, edge_count * sizeof(tetra));
  memcpy(t_arr + (count + edge_count), prev.t_arr, prev.len * sizeof(tetra));
  
  result.t_arr = t_arr;
  result.len = (count + edge_count + prev.len);
  free(cube_pts);
  free(prev_edges);
  free(prev.t_arr);
  
  printf("Finished dimension: %d,%d,%d\n",dim[0],dim[1],dim[2]);
  printf("  Previous dimension: %d,%d,%d\n",new_dim[0],new_dim[1],new_dim[2]);
  printf("  Tetraeders added: %zu\n", (count + edge_count));
  printf("  Tetraeders prev:  %zu\n", (prev.len));
  printf("  Amount of tetra: %zu\n", combination(len,4));
  printf("  Fraction sharp: %f\n", (float) result.len / combination(len,4));  
  printf("  Total tetraeders: %zu\n\n", result.len);
  return result;  
}
