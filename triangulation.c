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
#include "triangulation.h"

/*
 * Contains loads of functions which returns whether two objects
 * are disjoint (do not collide). Stuff for tet - tet tests and triangle - tet tests.
 */

void triangulation_free(ptriangulation triang) {
  free(triang->bound_tri);
  free(triang->tetra);
  free(triang);
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
 * Returns the relationship between two sets of vertices and this separation axis. 
 * Assuming we have simplices, we test if they touch face-to-face.
 * - Disjoint
 * - Touch
 * - Intersect
 */
int vert_vert_separation_axis(arr3 * vert1, int len1, arr3 * vert2, int len2,arr3 axis) {
  int min_1 = INT_MAX, max_1 = INT_MIN;
  int min_2 = INT_MAX, max_2 = INT_MIN;
  //Calculate projections of the first vertices set
  for (int i = 0; i < len1; i++) {
    int dot1 = dotArr3(axis, vert1[i]);
    if (dot1 < min_1) min_1 = dot1;
    if (dot1 > max_1) max_1 = dot1;
  }
  for (int i = 0; i < len2; i++) {
    int dot2 = dotArr3(axis, vert2[i]);
    if (dot2 < min_2) min_2 = dot2;
    if (dot2 > max_2) max_2 = dot2;
  }
  /*
   * Disjoint if either range [min_1, max_1] is before [min_2, max_2] or
   * [min_2, max_2] is before [min_1, max_1]
   */
  if (max_1 < min_2 || max_2 < min_1)
    return DISJOINT;
  else if (max_1 > min_2 && max_2 > min_1) //!(max_1 <= min_2 || max_2 <= min_1)
    return INTERSECT;
  else { 
    /*
     * The objects touch, this implies that some vertices of vert1 and vert2 lie in the sepation plane, while
     * the other vertices are separated by this plane. We have a different knd of touches:
     * face-to-face, edge-to-edge, vertex-to-vertex, and all mixed combinations. 
     * Note that the equation of a plane is given by dot(normal, vertex) = dot(normal, point_in_plane)
     * We now have that axis defines a (semi)-separating plane: vertices are split by the plane, or lie on the plane.
     * Now in our case dot(normal, point_in_plane) is given by either max_1 (if max_1 == min_2) or max_2 (if max_2 == min_1)
     *
     * Now we need to check if the vertices match. For example a face-to-face touch should have 3 points in the separating plane, with
     * all 3 matching. On the other hand, for a face-to-vertex touch, we should have 2 matching points.
     *
     * To do this, we calculate the amount of points in the plane for both simplices and the amount of points matched.
     * Next we calculate whether the amount of points matched equals the min of points in the plane of both simplices.
     */
    int dot_plane = (max_1 == min_2)? max_1 : max_2;
    int cnt_1 = vert_in_plane_count(axis, dot_plane, vert1, len1);
    int cnt_2 = vert_in_plane_count(axis, dot_plane, vert2, len2);
    int min_cnt = (cnt_1 < cnt_2)? cnt_1 : cnt_2;
    if (min_cnt == vert_vert_share_count(vert1,len1,vert2,len2))
      return DISJOINT;
    else
      return INTERSECT;//TOUCH
  }
}

/*
 * Applies the above function for a series of axis. Note that two objects can be only one of the three
 * states (DISJOINT/TOUCHING/INTERSECTING). Also, if for one axis we are TOUCHING, then the objects touch.
 * Similiary DISJOINT means the objects are disjoint. If we find no axis where we are TOUCHING/DISJOINT, then
 * we must INTERSECT
 */
int vert_vert_separation_axes(arr3 * vert1, int len1, arr3 * vert2, int len2,arr3 * axis,int len_axis) {
  for (int i = 0; i < len_axis; i++) {
    int sep = vert_vert_separation_axis(vert1, len1, vert2, len2, axis[i]);
    if (sep != INTERSECT) //We are TOUCHING/DISJOINT
      return sep;
  }
  return INTERSECT;
}

/*
 * Wrapper for the above functions in case we are checking tetra vs tetra or tri vs tetra
 */
int tet_tet_separation_axis(ptetra t1, ptetra t2, arr3 axis) {
  return vert_vert_separation_axis(t1->vertices, 4, t2->vertices, 4, axis);
}
int tet_tet_separation_axes(ptetra t1, ptetra t2, arr3 * axis, int len_axis) {
  return vert_vert_separation_axes(t1->vertices,4, t2->vertices,4, axis, len_axis);
}

int tri_tet_separation_axis(ptriangle t1, ptetra t2, arr3 axis) {
  return vert_vert_separation_axis(t1->vertices,3, t2->vertices,4, axis);
}
int tri_tet_separation_axes(ptriangle t1, ptetra t2, arr3 * axis, int len_axis) {
  return vert_vert_separation_axes(t1->vertices,3, t2->vertices,4, axis, len_axis);
}

/*
 * Returns whether this two tetrahedra are disjoint (in a face2face way). This means
 * that if two tetrahedra have a touching face they are disjoint iff they share this face.
 * If two tetrahedra have a touching edge they are disjoint iff they share this edge.
 */
int tet_tet_disjoint(ptetra t1, ptetra t2) {
  arr3 normals[4];
  /*
   * First check facet normals for separating axes.
   * This covers face-face and face-edge.
   */
  tetra_normals(t1, normals); 
  if (tet_tet_separation_axes(t1,t2,normals, 4) == DISJOINT)
    return DISJOINT;

  tetra_normals(t2, normals);
  if (tet_tet_separation_axes(t1,t2,normals, 4) == DISJOINT)
    return DISJOINT;

  /*
   * Check all cross produts of edges of both tetrahedra as separting axes.
   * This covers edge-edge.
   * 
   * The six edges of a tet are given by edges.
   * Indication how edge 'x' is made up.. for example edge[0] = vertex[1] - vertex[0]
   */
  int edges[6][2] = {{1,0},{2,0},{2,1}, //Base triangle
    {3,0},{3,1},{3,2}}; //Edges from base to apex
  //Test all the edges
  arr3 edge_a, edge_b;
  for (int i = 0; i < 6; i++) {
    subArr3(t1->vertices[edges[i][0]], t1->vertices[edges[i][1]], edge_a);
    for (int j = 0; j < 6; j++)
    {
      subArr3(t2->vertices[edges[j][0]], t2->vertices[edges[j][1]], edge_b);
      //Check crossproduct of edge_a and edge_b as seperating axis
      arr3 perp_axis;
      crossArr3(edge_a,edge_b, perp_axis);
      if (zeroArr3(perp_axis)) { 
        /*
         * Edges are parallel thus lie in some plane.
         * Calculate another vector in this plane and try the cross product
         * between those vectors as separating axis
         */
        subArr3(t1->vertices[edges[i][0]], t2->vertices[edges[j][0]], edge_b);
        crossArr3(edge_a,edge_b, perp_axis);
        if (zeroArr3(perp_axis)) //Still zero, all points lie on the same line, not a separation axis, IS TRUE?
          continue;
      }     
      if (tet_tet_separation_axis(t1,t2,perp_axis) == DISJOINT)
	return DISJOINT;
    }
  }
  return INTERSECT; //No separation axis found, thus we are not disjoint
}

int tri_tet_disjoint(ptriangle tri, ptetra tet) {
  arr3 normals[4];
  /*
   * First check facet normals for separating axes.
   * This covers face-face and face-edge.
   */
  triangle_normal(tri,normals[0]);
  if (tri_tet_separation_axis(tri,tet,normals[0]) == DISJOINT)
    return DISJOINT;

  tetra_normals(tet, normals); 
  if (tri_tet_separation_axes(tri,tet,normals, 4) == DISJOINT)
    return DISJOINT;

  /*
   * Check all cross produts of edges between triangle and tetrahedron as separating axis.
   * This covers edge-edge.
   * 
   * The six edges of a tet are given by edges.
   * Indication how edge 'x' is made up.. for example edge[0] = vertex[1] - vertex[0]
   */
  int edges[6][2] = {{1,0},{2,0},{2,1}, //Base triangle (first three)
    {3,0},{3,1},{3,2}}; //Edges from base to apex
  //Test all the edges
  arr3 edge_a, edge_b;

  for (int i = 0; i < 3; i++) {
    subArr3(tri->vertices[edges[i][0]], tri->vertices[edges[i][1]], edge_a);
    for (int j = 0; j < 6; j++)
    {
      subArr3(tet->vertices[edges[j][0]], tet->vertices[edges[j][1]], edge_b);
      //Check crossproduct of edge_a and edge_b as seperating axis
      arr3 perp_axis;
      crossArr3(edge_a,edge_b, perp_axis);
      if (zeroArr3(perp_axis)) { 
        /*
         * Edges are parallel thus lie in some plane.
         * Calculate another vector in this plane and try the cross product
         * between those vectors as separating axis
         */
        subArr3(tri->vertices[edges[i][0]], tet->vertices[edges[j][0]], edge_b);
        crossArr3(edge_a,edge_b, perp_axis);
        if (zeroArr3(perp_axis)) //Still zero, all points lie on the same line, not a separation axis, IS TRUE?
          continue;
      }     
      if (tri_tet_separation_axis(tri,tet, perp_axis) == DISJOINT)
	return DISJOINT;
    }
  }
  return INTERSECT; //No separation axis found, thus we are not disjoint
}


int tet_tet_list_disjoint(ptetra tet, ptetra  tet_list, int list_len) {
  for (int i = 0; i < list_len; i++) 
    if (tet_tet_disjoint(tet,tet_list+ i) == INTERSECT)
      return INTERSECT;  

  //No intersection, so tetrahedron must be disjoint from the list
  return DISJOINT;
}



int tet_triangulation_disjoint(ptetra tet, ptriangulation triang) {
  for (size_t i = 0; i < triang->tetra_len; i++)
    if (tet_tet_disjoint(tet, triang->tetra + i) == INTERSECT)
      return INTERSECT;

  return DISJOINT;
}

int tri_triangulation_disjoint(ptriangle tri, ptriangulation triang) {
  for (size_t i = 0; i < triang->tetra_len; i++)
    if (tri_tet_disjoint(tri, triang->tetra + i) == INTERSECT)
      return INTERSECT;

  return DISJOINT;
}

/*
 * Returns -1 if not found
 */
int tri_arr_idx(ptriangle tri, ptriangle list, size_t list_len) {
  for (size_t i = 0; i < list_len; i++)
    if (tri_tri_equal(tri, list + i))
      return i;
  return -1;
}


void rem_boundary_triangulation(int rem_bound, ptriangulation result) {
  result->bound_len--;
  ptriangle new_bound = malloc(result->bound_len * sizeof(triangle));
  memcpy(new_bound,result->bound_tri,rem_bound * sizeof(triangle));
  memcpy(new_bound + rem_bound, result->bound_tri + rem_bound + 1, (result->bound_len - rem_bound) * sizeof(triangle));
  free(result->bound_tri);
  result->bound_tri = new_bound;
}

int update_boundary_triangulation(arr3 v1, arr3 v2, arr3 v3, ptriangulation triang) {
  triangle tri = (triangle) {{{v1[0],v1[1],v1[2]}, {v2[0],v2[1],v2[2]}, {v3[0],v3[1],v3[2]}}};
  int idx = (tri_arr_idx(&tri, triang->bound_tri, triang->bound_len));
  if (idx > -1) //Wanting to add a boundary again, actually remove!
  {
    rem_boundary_triangulation(idx, triang);
    return -2;
  }
  if (triangle_boundary_cube(&tri,triang->dim)) //Triangle on the boundary of the cube
    return -1;

  printf("New boundary for triangulation:\n");
  print_triangle(&tri);
  triang->bound_len++;
  triang->bound_tri = realloc(triang->bound_tri,triang->bound_len * sizeof(triangle));
  triang->bound_tri[triang->bound_len - 1] = tri;
  return 1;
}



void add_tet_triangulation(ptetra tet, ptriangulation result) {
  result->tetra_len++;
  result->tetra = realloc(result->tetra,result->tetra_len * sizeof(tetra));
  result->tetra[result->tetra_len - 1] = *tet;
  //Tetrahedron has 4 facets, check if facets are on boundary of triang / cube.
  //0,1,2
  update_boundary_triangulation(tet->vertices[0], tet->vertices[1], tet->vertices[2], result);
  //1,2,3
  update_boundary_triangulation(tet->vertices[1], tet->vertices[2], tet->vertices[3], result);
  //0,2,3
  update_boundary_triangulation(tet->vertices[0], tet->vertices[2], tet->vertices[3], result);
  //0,1,3
  update_boundary_triangulation(tet->vertices[0], tet->vertices[1], tet->vertices[3], result);
} 


void filter_tet_list_disjoint_triangulation(ptetra   list, size_t * list_len, ptriangulation triang) {
  size_t c = 0;
  for (size_t i = 0; i< *list_len; i++)
    if (tet_triangulation_disjoint(list+ i, triang)) {
      list[c] = list[i];
      c++;
    }

  *list_len = c;
}

int consistent_triangulation(ptriangulation triang, facet_acute_data * data) {
  int consistent = 1;
  for (size_t i = 0; i < triang->bound_len; i++) {
    if (!facet_conform(triang->bound_tri + i, data))  {
      printf("Somehow we have a non conform facet on the boundary.. WTF?\n");
      print_triangle(triang->bound_tri + i);
      consistent = 0;
    }
    if (!tri_list_contains(&data->data->list, triang->bound_tri + i)) {
      printf("Somehow the bondary triangle is not in the list anymore.. WTF?\n");
      print_triangle(triang->bound_tri + i);
      consistent = 0;
    }
  }
  return consistent;
}
size_t filter_tri_list_disjoint_tet(ptriangulation triang,facet_acute_data * parameters, tri_list * list, ptetra tet) {
  size_t cnt = tri_list_count(list);
  triangle cur_tri;
  size_t dim_size = (list->dim + 1) * (list->dim + 1) * (list->dim + 1);
  size_t i,j,k;
  int l;

  #pragma omp parallel for schedule(dynamic,list->dim) private(j,k,i,l,cur_tri) 
  for (i = 0; i < dim_size; i++) {
    vertex_from_index_cube(i,list->dim, cur_tri.vertices[0]);

    for (j = i; j < dim_size; j++) {
      vertex_from_index_cube(j, list->dim, cur_tri.vertices[1]);

      for (l = list->t_arr[i][j-i].len - 1; l >= 0; l--) {  //Loop over all triangles (i,j,*)
        k = list->t_arr[i][j-i].p_arr[l] +  j;

        vertex_from_index_cube(k, list->dim, cur_tri.vertices[2]);
        //Cur_tri is the l-th triangle with points(i,j,*)
        if (!tri_tet_disjoint(&cur_tri, tet)) {//If not disjoint with new tetrahedron, delete
          tri_list_remove(list, &cur_tri); //Thread safe, as only one processor acces [i][j]
	  /*
          if (!consistent_triangulation(triang,parameters)) {
            printf("Somehoew removing this triangle resulted in incorrect triangulation\n");
            print_triangle(&cur_tri);
            exit(0);
          }
	  */
        }
      }
    }
  }
  return (cnt - tri_list_count(list));
}

ptriangulation triangulate_cube(tri_list * list) {

  ptriangulation result = calloc(sizeof(triangulation), 1);
  result->dim = list->dim;

  /* Try to find a start facet on the z = 0 plane */
  ptriangle  start_facet = calloc(sizeof(triangle), 1);
  int found_start = 0;
  while (!found_start) { 
    /*
     * Random triangle with vertices:
     * (0,0,0), (?,0,0), (?,?,0)
     */
    start_facet->vertices[1][0] = rand() %list->dim;
    start_facet->vertices[2][0] = rand() %list->dim;
    start_facet->vertices[2][1] = rand() %list->dim;

    if (!triangle_acute(start_facet))
      continue;
    found_start =  (tri_list_contains(list, start_facet));
  }

  facet_acute_data parameters;
  data_list data;
  data.mode = DATA_TRI_LIST;
  data.list = *list;

  cube_points cube = gen_cube_points(list->dim);
  parameters.cube = &cube;
  parameters.boundary_func = &triangle_boundary_cube;
  parameters.data = &data;
  parameters.store_tetra = 1;

  //Start triangle (0,0,0), (rand,0,0), (rand,rand,0)
  printf("Starting triangulation with facet:\n");
  print_triangle(start_facet);
  printf("Triangle acute? %d\n", triangle_acute(start_facet));
  result->bound_len = 1;
  result->bound_tri = start_facet;

  //While we have triangles on the boundary..
  while (result->bound_len > 0) {
    /*
     * We are going to add a tetrahedron on the boundary triangle.
     * To do so, we select a random triangle on the boundary. Then we generate all the
     * acute tetrahedra (above and below) with facets in our possible list.
     * From this list we remove all the tetrahedrons that intersect with our current triangulation.
     * Then we add a random tetrahedron to our triangulation, update the conform list and repeat.
     */
    int rand_bound = rand() % result->bound_len;
    printf("\n\nTotal amount of triangles left:%zu\nExpanding triangulation at boundary triangle: \n", tri_list_count(list));
    print_triangle(result->bound_tri + rand_bound);

    //Calculate the conform tetrahedrons above and below
    facet_conform(&result->bound_tri[rand_bound], &parameters);

    // Concentanate ... Do this in face_cube_acute already I guess 
    size_t list_len = parameters.tet_above_len + parameters.tet_below_len;
    ptetra tet_list = malloc(list_len * sizeof(tetra));
    memcpy(tet_list                             , parameters.tet_above, parameters.tet_above_len * sizeof(tetra));
    memcpy(tet_list + parameters.tet_above_len, parameters.tet_below, parameters.tet_below_len * sizeof(tetra));
    free(parameters.tet_below); free(parameters.tet_above);

    printf("Total amount of tetrahedrons found: %zu\n", list_len);
    //Remove all the tetrahedrons that collide with current triangulation.
    filter_tet_list_disjoint_triangulation(tet_list, &list_len,result);

    printf("Amount of tetrahedrons left after filtering: %zu\n\n",list_len);
    if (list_len == 0) {
      printf("Waarom is deze lijst nu al fucking leeggefilterd?\n");
      exit(0);
      printf("Dead end, helaas pindakaas. Got to %zu\n", result->tetra_len);
      free(cube.points);
      triangulation_free(result);
      return NULL;
    }


    int rand_tet = rand() % list_len;
    /*
     * Add one of the tetrahedrons to the triangulation.
     * This removes all the boundary triangles that are covered by this tetrahedron
     */
    printf("Adding the following tetra to the triangulation\n");
    print_tetra(tet_list + rand_tet);
    printf("\n\n");
    add_tet_triangulation(tet_list + rand_tet, result);
    printf("\nNew amount of tet in triang %zu. New amount of boundaries %zu.\n\n", result->tetra_len, result->bound_len);

    //Consistency check
    consistent_triangulation(result, &parameters);

    /*
     * Remove all the triangles in list that collide with this tetrahedron.
     * Maybe save a list of triangles we have removed? This allows us to restore them,
     * in case we get a dead end.
     */
    printf("Removing all triangles that intersect with the new tetra\n");
    size_t removed =  filter_tri_list_disjoint_tet(result,&parameters,list, tet_list+rand_tet);
    printf("Removed %zu triangles not disjoint with new tetrahedron\n\n", removed);

    consistent_triangulation(result, &parameters);
    /*
     * Conform the resulting set of triangles. We probably do not want to do this step every time,
     * what is a good measure?
     */
    
    //facets_conform(&data, NULL);

    free(tet_list);

  }
  free(cube.points);
  printf("Triangulation has length of %zu\n", result->tetra_len);
  return result;
}

/*
 * Semi-random triangulation of the cube. Start with the fundamental index list.
 */
//ptriangulation triangulate_cube_random(int dim, tri_index_list   * fund_ind_list) {
 // ptriangulation result = calloc(sizeof(triangulation), 1);
 // /*
 // result->dim = dim;

 // /* The total index list */
 // tri_index_list ind_list;
 // ind_list.len = fund_ind_list->len * 48
 //   tri_index * t_arr = malloc(sizeof(tri_index) * ind_list.len);

 // /*
 //  * Apply symmetry on each triangle from the fundamental list, to get the
 //  * entire list
 //  */
 // for (size_t i = 0; i < fund_ind_list.len; i++) {



 // }
 // /* Try to find a start facet on the z = 0 plane */
 // ptriangle      start_facet = calloc(sizeof(triangle), 1);
 // int found_start = 0;
 // while (!found_start) { 
 //   /*
 //    * Random triangle with vertices:
 //    * (0,0,0), (?,0,0), (?,?,0)
 //    */
 //   start_facet->vertices[1][0] = rand() % dim;
 //   start_facet->vertices[2][0] = rand() % dim;
 //   start_facet->vertices[2][1] = rand() % dim;
 //   facet_conform
 //     found_start = facet_cube_acute(start_facet, &parameters, FACET_ACUTE);
 // }



 // facet_acute_data parameters;
 // cube_points cube = gen_cube_points(dim);
 // parameters.cube = &cube;

 // //Start triangle (0,0,0), (rand,0,0), (rand,rand,0)
 // printf("Found acute facet:\n");
 // print_triangle(start_facet);
 // add_boundary_triangulation(start_facet, result);
 // free(start_facet);
 // while (result->bound_len > 0) {
 //   int rand_bound = rand() % result->bound_len;
 //   facet_cube_acute(&result->bound_tri[rand_bound], &parameters, FACET_ACUTE_TETRA); 

 //   // Concentanate ... Do this in face_cube_acute already I guess 
 //   size_t list_len = parameters.tetra_above_len + parameters.tetra_below_len;
 //   ptetra tet_list = malloc(list_len * sizeof(tetra));
 //   memcpy(tet_list                             , parameters.tetra_above, parameters.tetra_above_len * sizeof(tetra));
 //   memcpy(tet_list + parameters.tetra_above_len, parameters.tetra_below, parameters.tetra_below_len * sizeof(tetra));
 //   free(parameters.tetra_below); free(parameters.tetra_above);
 //   // Filtering.. Do this in facet_cube_acute already? 
 //   filter_tetra_list_disjoint(&tet_list, &list_len,result);
 //   if (list_len == 0) {
 //     printf("Waarom is deze lijst nu al fucking leeggefilterd?\n");
 //     exit(0);
 //     printf("Dead end, helaas pindakaas. Got to %zu\n", result->tetra_len);
 //     free(cube.points);
 //     triangulation_free(result);
 //     return NULL;
 //   }

 //   int rand_tet = rand() % list_len;
 //   add_tet_triangulation(tet_list + rand_tet, result); //Add this tetrahedron to the triangulation
 //   rem_boundary_triangulation(rand_bound,result);
 // }
 // free(cube.points);
 // printf("Triangulation has length of %d\n", result->tetra_len);
 // */
//  return result;
//}

  /*
int tri_facet_boundary_triangulation(arr3 v1, arr3 v2, arr3 v3, ptriangulation result, ptriangle search_tri) {
   *search_tri = (triangle) {{{v1[0],v1[1],v1[2]}, {v2[0],v2[1],v2[2]}, {v3[0],v3[1],v3[2]}}};
   if (triangle_boundary_cube(search_tri, result->dim))
   return 1;
   tri_index search_facet;
   triangle_to_index_cube((*search_tri), result->dim_mult, search_facet);
   for (size_t i = 0; i < result->bound_len; i++) {
   tri_index bound_facet;
   triangle_to_index_cube(result->bound_tri[i],result->dim_mult, bound_facet);
   if (bound_facet[0] == search_facet[0] && bound_facet[1] == search_facet[1] &&
   bound_facet[2] == search_facet[2])
   return 1;
   }
   return 0;
  return 0;
}
   */
