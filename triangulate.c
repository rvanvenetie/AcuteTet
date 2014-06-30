#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "limits.h"
#include "triangle.h"
#include "tetraeder.h"
#include "triangulate.h"

void triangulation_free(ptriangulation triang) {
  free(triang->boundaries);
  free(triang->tetraeders);
}

void tetra_edges(ptetra tet, arr3 * edges) {
  //Edges of base triangle
  subArr3(tet->vertices[1], tet->vertices[0], edges[0]);
  subArr3(tet->vertices[2], tet->vertices[0], edges[1]);
  subArr3(tet->vertices[2], tet->vertices[1], edges[2]);
  //Now only need the three edges from base to apex 
  subArr3(tet->vertices[3], tet->vertices[0], edges[3]);
  subArr3(tet->vertices[3], tet->vertices[1], edges[4]); 
  subArr3(tet->vertices[3], tet->vertices[2], edges[5]);     
}
/*
 * Test for non-strict overlapping
 */
int tetra_disjoint_axis(ptetra t1, ptetra t2, arr3 axis) {
  int min_1 = INT_MAX, max_1 = INT_MIN;
  int min_2 = INT_MAX, max_2 = INT_MIN;
  for (int i = 0; i < 4; i++) {
    int dot1 = dotArr3(axis, t1->vertices[i]);
    if (dot1 < min_1) min_1 = dot1;
    if (dot1 > max_1) max_1 = dot1;

    int dot2 = dotArr3(axis, t2->vertices[i]);
    if (dot2 < min_2) min_2 = dot2;
    if (dot2 > max_2) max_2 = dot2;
  }
  //Now we know the scalar projections on the axis, test for overlap
  
  /*
   * Disjoint if range A is completely after range B, or that range A is completely before B.
   */
  return (min_1 >= max_2 || max_1 <= min_1);
}
/*
 * If disjoint, we return true, as this means we have found a separating plane
 */
int tetra_disjoint_axes(ptetra t1, ptetra t2, arr3 * axes, int len) {
  for (int i = 0; i < len; i++) {
    if (tetra_disjoint_axis(t1,t2,axes[i]))
      return 1;
  }
  return 0;
}
int tetra_tetra_disjoint(ptetra t1, ptetra t2) {
  arr3 normals[4];
  tetra_normals(t1, normals); //First check facets for separating axis
  if (tetra_disjoint_axes(t1,t2,normals, 4))
    return 1;
  tetra_normals(t2, normals);
  if (tetra_disjoint_axes(t1,t2,normals,4))
    return 1;
  
  //Indication how edge i is made up.. for example edge[0] = vertex[1] - vertex[0]
  int edges[6][2] = {{1,0},{2,0},{2,1}, //Base triangle
                   {3,0},{3,1},{3,2}}; //Edges from base to apex
  //Test all the edges
  arr3 edge_a, edge_b;
  for (int i = 0; i < 6; i++) {
    subArr3(t1->vertices[edges[i][0]], t1->vertices[edges[i][1]], edge_a);
    for (int j = 0; j < 6; j++)
    {
      subArr3(t2->vertices[edges[j][0]], t2->vertices[edges[j][1]], edge_b);
      //Find normal to plane parallel to edge_a and edge_b
      arr3 perp_axis;
      crossArr3(edge_a,edge_b, perp_axis);
      if (zeroArr3(perp_axis)) { //Edges lie in the same plane, try vector from a to b
        subArr3(t1->vertices[edges[i][0]], t2->vertices[edges[j][0]], edge_b);
        crossArr3(edge_a,edge_b, perp_axis);
        if (zeroArr3(perp_axis)) //Still zero, all points lie on the same line, not a separation axis, IS TRUE?
          continue;
      }     
      if (tetra_disjoint_axis(t1,t2, perp_axis))
        return 1;
    }
  }
  return 0; //No separation axis found, thus we are not disjoint
}

int tetra_triangulation_disjoint(ptetra tet, ptriangulation triang) {
  for (size_t i = 0; i < triang->tetra_len; i++)
    if (!tetra_tetra_disjoint(tet, &triang->tetraeders[i])) {
    
    /*
      print_tetra(tet);
      printf("Not disjoint with: \n");
      print_tetra(triang->tetraeders + i);
      */
      return 0;
    }
  return 1;
}

void filter_tetra_list_disjoint(ptetra *  list, size_t * list_len, ptriangulation triang) {
  /*
  size_t c = 0;
  for (size_t i = 0; i< *list_len; i++)
    if (tetra_triangulation_disjoint((*list) + i, triang)) {
      (*list)[c] = (*list)[i];
      c++;
    }

  *list_len = c;
  *list = realloc(*list, c * sizeof(tetra));  
  */
}

void add_boundary_triangulation(ptriangle triang, ptriangulation result) {
  /*
  printf("New boundary for triangulation:\n");
  print_triangle(triang);
  result->bound_len++;
  result->boundaries = realloc(result->boundaries,result->bound_len * sizeof(triangle));
  result->boundaries[result->bound_len - 1] = *triang;
  */
}

void rem_boundary_triangulation(int rem_bound, ptriangulation result) {
  /*
  result->bound_len--;
  ptriangle new_bound = malloc(result->bound_len * sizeof(triangle));
  memcpy(new_bound,result->boundaries,rem_bound * sizeof(triangle));
  memcpy(new_bound + rem_bound, result->boundaries + rem_bound + 1, (result->bound_len - rem_bound) * sizeof(triangle));
  free(result->boundaries);
  result->boundaries = new_bound;
  */
}

int facet_boundary_triangulation(arr3 v1, arr3 v2, arr3 v3, ptriangulation result, ptriangle search_tri) {
  /*
  *search_tri = (triangle) {{{v1[0],v1[1],v1[2]}, {v2[0],v2[1],v2[2]}, {v3[0],v3[1],v3[2]}}};
  if (triangle_boundary_cube(search_tri, result->dim))
    return 1;
  tri_index search_facet;
  triangle_to_index_cube((*search_tri), result->dim_mult, search_facet);
  for (size_t i = 0; i < result->bound_len; i++) {
    tri_index bound_facet;
    triangle_to_index_cube(result->boundaries[i],result->dim_mult, bound_facet);
    if (bound_facet[0] == search_facet[0] && bound_facet[1] == search_facet[1] &&
        bound_facet[2] == search_facet[2])
      return 1;
  }
  return 0;
  */
}

void add_tet_triangulation(ptetra tet, ptriangulation result) {
  /*
  result->tetra_len++;
  result->tetraeders = realloc(result->tetraeders,result->tetra_len * sizeof(tetra));
  result->tetraeders[result->tetra_len - 1] = *tet;
  printf("Maar drie prints mogen hier komen");
  //Tetrahedron has 4 facets, check if facets are on boundary of triang / cube.
  triangle facet;
  //0,1,2
  if (!facet_boundary_triangulation(tet->vertices[0], tet->vertices[1], tet->vertices[2], result, &facet)) //Not on boundary
    add_boundary_triangulation(&facet,result);
  //1,2,3
  if (!facet_boundary_triangulation(tet->vertices[1], tet->vertices[2], tet->vertices[3], result, &facet)) //Not on boundary
    add_boundary_triangulation(&facet,result);
  //0,2,3
  if (!facet_boundary_triangulation(tet->vertices[0], tet->vertices[2], tet->vertices[3], result, &facet)) //Not on boundary
    add_boundary_triangulation(&facet,result);
  //0,1,3
  if (!facet_boundary_triangulation(tet->vertices[0], tet->vertices[1], tet->vertices[3], result, &facet)) //Not on boundary
    add_boundary_triangulation(&facet,result);
  */
} 

ptriangulation triangulate_cube_random(arr3 dim) {
  /*
  ptriangulation result = calloc(sizeof(triangulation), 1);
  result->dim[0] = dim[0];
  result->dim[1] = dim[1];
  result->dim[2] = dim[2];
  
  result->dim_mult[0] = (dim[0] + 1) * (dim[1] + 1);
  result->dim_mult[1] = (dim[0] + 1);
  
  facet_acute_data parameters;
  cube_points cube = gen_cube_points(dim);
  parameters.cube = &cube;
  
  //Start triangle (0,0,0), (rand,0,0), (rand,rand,0)
  ptriangle      start_facet = calloc(sizeof(triangulation), 1);
  int start_acute = 0;
  while (!start_acute) { //Find start facet on z=0 plane
    start_facet->vertices[1][0] = rand() % dim[0];
    start_facet->vertices[2][0] = rand() % dim[0];
    start_facet->vertices[2][1] = rand() % dim[1];
    start_acute = facet_cube_acute(start_facet, &parameters, FACET_ACUTE);
  }
  printf("Found acute facet:\n");
  print_triangle(start_facet);
  add_boundary_triangulation(start_facet, result);
  free(start_facet);
  while (result->bound_len > 0) {
    int rand_bound = rand() % result->bound_len;
    facet_cube_acute(&result->boundaries[rand_bound], &parameters, FACET_ACUTE_TETRA); 
    
    // Concentanate ... Do this in face_cube_acute already I guess 
    size_t list_len = parameters.tetra_above_len + parameters.tetra_below_len;
    ptetra tet_list = malloc(list_len * sizeof(tetra));
    memcpy(tet_list                             , parameters.tetra_above, parameters.tetra_above_len * sizeof(tetra));
    memcpy(tet_list + parameters.tetra_above_len, parameters.tetra_below, parameters.tetra_below_len * sizeof(tetra));
    free(parameters.tetra_below); free(parameters.tetra_above);
    // Filtering.. Do this in facet_cube_acute already? 
    filter_tetra_list_disjoint(&tet_list, &list_len,result);
    if (list_len == 0) {
      printf("Waarom is deze lijst nu al fucking leeggefilterd?\n");
      exit(0);
      printf("Dead end, helaas pindakaas. Got to %zu\n", result->tetra_len);
      free(cube.points);
      triangulation_free(result);
      return NULL;
    }
    
    int rand_tet = rand() % list_len;
    add_tet_triangulation(tet_list + rand_tet, result); //Add this tetrahedron to the triangulation
    rem_boundary_triangulation(rand_bound,result);
  }
  free(cube.points);
  printf("Triangulation has length of %d\n", result->tetra_len);
  return result;
  */
}
