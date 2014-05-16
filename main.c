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


#define USE_FILE 1
#define FACE2FACE 1
#define LOOP


int facet_tetra_list(ptriangle triang, arr3 new_vertex, tri_mem_list * acute_list) {
  return (mem_list_get_sym_fund(acute_list, triang->vertices[0], triang->vertices[1], new_vertex) &&
          mem_list_get_sym_fund(acute_list, triang->vertices[0], triang->vertices[2], new_vertex) &&
          mem_list_get_sym_fund(acute_list, triang->vertices[1], triang->vertices[2], new_vertex));
}


int facet_cube_acute(ptriangle triang, cube_points * cube, tri_mem_list * acute_list) {
  //Check if triangle is acute..
  if (!mat3_triangle_acute(triang->vertices)) //Can replace with code below
    return 0;
  mat3 P;
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
  
  
  int bound_triang = triangle_boundary(triang,cube->dim); //Boundary plane only needs acute tetra on 1 side
  int dotprod;
  int found_above = 0, found_below = 0;
  
  for (size_t i = 0; i < cube->len; i++){
    dotprod = dotArr3(cube->points[i], tri_normal); //Calculate angle between normal and vector from plane to cube_pt
    //Check if current point lies on side of the triangle which isn't sharp yet and check if point lies in the
    //prism of possible sharp tetrahedron
    
    if (((dotprod > tri_d && !found_above) ||  //Dotprod > tri_d implies that the point lies "above" the triangle (same side as normal) 
         (dotprod < tri_d && !found_below))     //Dotprod < tri_d implies lies below
          && //Check if point lies in the prism
          dotArr3(cube->points[i],side_normals[0]) < side_d[0] &&
          dotArr3(cube->points[i],side_normals[1]) < side_d[1] &&
          dotArr3(cube->points[i],side_normals[2]) < side_d[2]) 
      {
      if (acute_list && !facet_tetra_list(triang, cube->points[i], acute_list))
        continue;  
      tetra test_tetra;
      memcpy(test_tetra.vertices, &cube->points[i], sizeof(arr3));
      memcpy(test_tetra.vertices + 1, triang->vertices, 3 * sizeof(arr3)); 
      if (tetra_acute(&test_tetra)) {
        if (dotprod > tri_d) 
          found_above = 1;
        else 
          found_below = 1;
        if ((found_above && found_below) || bound_triang)
          return 1;
      }
    }
  }
  return 0;  
}


tri_mem_list facets_cube_acute(int dim) {
  tri_mem_list result = mem_list_init_fund(dim);
  arr3 tmpdim = {dim,dim,dim};
  cube_points fund = gen_fund_points(dim);
  cube_points cube = gen_cube_points(tmpdim);
  size_t comb_len;
  Dindex * comb = combinations_list(cube.len, 2, &comb_len);
  
  size_t j,k,i,l;
  triangle cur_tri;
  
  #pragma omp parallel for schedule(guided) private(l,j,k,i,cur_tri)

  for (l = 0; l < comb_len; l++) {//Loop through all combinations
    j = comb[l*2];
    k = comb[l*2+1];
    for (i = 0; i < fund.len; i++){    
      cur_tri = (triangle) {{{fund.points[i][0],fund.points[i][1],fund.points[i][2]},
                             {cube.points[j][0],cube.points[j][1],cube.points[j][2]},
                             {cube.points[k][0],cube.points[k][1],cube.points[k][2]}}};
      if (facet_cube_acute(&cur_tri,&cube,NULL))
        mem_list_set_sym_fund(&result, &cur_tri);
    }
  }
  free(comb);
  free(cube.points);
  free(fund.points);
  return result;
}


void facets_face2face(tri_mem_list * acute_list){
  cube_points fund = gen_fund_points(acute_list->dim[0]);
  cube_points cube = gen_cube_points(acute_list->dim);
  size_t comb_len;
  Dindex * comb = combinations_list(cube.len, 2, &comb_len);
  int changed = 1;
  size_t l,i,j,k;
  tri_index indices;
  unsigned short idx2, idx3;
  triangle cur_tri;
  while (changed) {
    changed = 0;
    printf("Starting face2face loop with %zu acute triangles, from thread %d\n"
          , mem_list_count(acute_list), omp_get_thread_num());
    #pragma omp parallel for schedule(guided) private(l,j,k,i,cur_tri, idx2,idx3,indices) 
    for (l = 0; l < comb_len; l++) {//Loop through all combinations
      j = comb[l*2];
      k = comb[l*2+1];
      idx2 = vertex_to_index(cube.points[j], acute_list->dim_mult);
      idx3 = vertex_to_index(cube.points[k], acute_list->dim_mult);
      for (i = 0; i < fund.len; i++){ //Against all fundamental points
        vertices_unique_fund(vertex_to_index(fund.points[i],acute_list->dim_mult), idx2, idx3, acute_list, indices); 
        if (!GMI(acute_list->t_arr,indices)) //Check if this index is still acute
          continue;
        cur_tri = (triangle) {{{fund.points[i][0],fund.points[i][1],fund.points[i][2]},
                             {cube.points[j][0],cube.points[j][1],cube.points[j][2]},
                             {cube.points[k][0],cube.points[k][1],cube.points[k][2]}}};
        if (!facet_cube_acute(&cur_tri,&cube,acute_list)) { //remove from list
          changed = 1;
          mem_list_clear_sym_fund(acute_list, &cur_tri);
        }
      }
    } 
  }
  
  
  free(comb); 
  free(cube.points);
  free(fund.points);
}

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

int main(void){
  test_multithread();
  //exit(00);
  char filename[70];
  //arr3 dim;// =  {7,7,7};
  tri_mem_list fund_list, face_list;
  //triangle_list tri_list;
  //tri_index_list idx_list;
  //tetra_list tet_list;
  
  double time_start,time_end;
  
 // tet_list = acute_tetrahedra_recur(dim);
 // exit(0);
  int i;
  #ifdef LOOP
  for (i = 1;i <21; i++) 
  #endif
  {
    //arr2 dim_mult = {(i+1)*(i+1), i+1};
    //dim[0] = i; dim[1] = i; dim[2] = i;
    printf("Dimension: %d. Gathering fundamental data.\n", i);
    sprintf(filename,"/var/scratch/rvveneti/fund_data_%d.tet",i);
    if (USE_FILE && mem_list_from_file(&fund_list, filename))
      printf("Succesfully data from file%s\n", filename);
    else {
      time_start = omp_get_wtime();
      fund_list = facets_cube_acute(i);
      time_end   = omp_get_wtime();
      printf("Calculated data in %f seconds\n", time_end - time_start);  
      if (!mem_list_to_file(&fund_list, filename)) 
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
      facets_face2face(&face_list);
      time_end   = omp_get_wtime();
      printf("Calculated data in %f seconds\n", time_end   - time_start);  
      printf("Memory before clean up: %zu\n", mem_list_memory(&face_list));
      mem_list_clean(&face_list);
      if (!mem_list_to_file(&face_list, filename))
        printf("Failed to save face2face data to file: %s\n", filename);
    }
    printf("Amount of sharp facets after face2face filter: %zu\n", mem_list_count(&face_list));    
    printf("Total memory used: %zu\n\n", mem_list_memory(&face_list));
    mem_list_free(&face_list);
  }
}
