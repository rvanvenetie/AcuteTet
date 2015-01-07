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

triangulation triangulation_init(int dim) {
  triangulation result = {0};
  result.dim = dim;
  return result;
}

void triangulation_free(ptriangulation triang) {
  free(triang->bound_tri);
  free(triang->tetra);
}

float triangulation_volume(ptriangulation triang) {
  int vol_tot = 0;
  for (size_t i = 0; i < triang->tetra_len; i++)
    vol_tot += tetra_volume(triang->tetra + i);

  return ((float) vol_tot / 6.0);
}

void triangulation_print(ptriangulation triang) {
  printf("Triangulation. Dimension = %d. Boundary triangles = %zu. Tetra = %zu. Volume covered %f.\n",\
      triang->dim,triang-> bound_len, triang->tetra_len, triangulation_volume(triang) / ((triang->dim+1)*(triang->dim+1)*(triang->dim+1)));
}

int triangulation_to_file(ptriangulation triang, char * filename) {
  FILE * stream;
  stream = fopen(filename, "wb");
  if (stream == NULL)
    return 0;
  //Write the struct to the file
  if (fwrite(triang, sizeof(triangulation), 1, stream) < 1)
    return 0;
  if (fwrite(triang->bound_tri, sizeof(triangle), triang->bound_len, stream) < triang->bound_len)
    return 0;
  if (fwrite(triang->tetra, sizeof(tetra), triang->tetra_len, stream) < triang->tetra_len)
    return 0;

  return 1;
}

int triangulation_from_file(ptriangulation triang, char * filename) {
  FILE * stream;
  stream = fopen(filename, "rb");
  if (stream == NULL)
    return 0;

  //Reads the entire struct from file
  if (fread(triang, sizeof(triangulation), 1, stream) < 1)
    return 0;
  triang->bound_tri = malloc(sizeof(triangle) * triang->bound_len);
  if (fread(triang->bound_tri, sizeof(triangle), triang->bound_len, stream) < triang->bound_len)
    return 0;
  triang->tetra = malloc(sizeof(tetra) * triang->tetra_len);
  if (fread(triang->tetra, sizeof(tetra), triang->tetra_len, stream) < triang->tetra_len)
    return 0;

  return 1;
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


void filter_tet_list_disjoint_triangulation(ptetra   list, unsigned short * list_len, ptriangulation triang) {
  unsigned short c = 0;
  for (unsigned short i = 0; i< *list_len; i++)
    if (tet_triangulation_disjoint(list+ i, triang)) {
      list[c] = list[i];
      c++;
    }

  *list_len = c;
}

void tet_list_min_max_volume(ptetra list, unsigned short list_len, unsigned short * max, unsigned short * min) {
  *min = 0;
  *max = 0;
  int vol_max = INT_MIN;
  int vol_min = INT_MAX;

  for (unsigned short i = 0; i < list_len; i++) 
  {
    int vol_cur = tetra_volume(list + i);
    if (vol_cur < vol_min) {
      vol_min = vol_cur;
      *min = i;
    }

    if (vol_cur > vol_max) {
      vol_max = vol_cur;
      *max = i;
    }
  }
}

int triangulation_consistent(ptriangulation triang, facet_acute_data * data) {
  int consistent = 1;
  int store_acute_ind = data->store_acute_ind;
  data->store_acute_ind = 0;
  for (size_t i = 0; i < triang->bound_len; i++) {
    if (!mem_list_cube_contains(&data->data->mem_list, triang->bound_tri + i)) // Boundary triangle removed from the data
    {
      printf("Triangle on the boundary got removed from the data.\n");
      consistent = 0;
      break;
    }
    if (!facet_conform(triang->bound_tri + i, data))  { //Triangle on the boundary not conform any longer
      printf("Triangle on the boundary not a conform facet anymore\n");
      consistent = 0;
      break;
    }
  }
  data->store_acute_ind = store_acute_ind;
  return consistent;
}
/*
 * Everything is given by indices.
 * The function below constructs an array of tetra given by the base_triangle and the apices.
 * It adds all the facets != base_triangle to a given list.
 *
 * 2D Array of locks given as parameter, in case this function gets called in OpenMP
 */
size_t facets_add_cnt;
#pragma omp threadprivate(facets_add_cnt)


void facets_tetra_list(tri_list * list, tri_index base_idx, vert_index * ind_apices, int apices_len, omp_lock_t ** locks) {
  for (int i = 0; i < apices_len; i++)
  {
    tri_index  facets[3] = {{base_idx[0], base_idx[1], ind_apices[i]},
      {base_idx[1], base_idx[2], ind_apices[i]},
      {base_idx[0], base_idx[2], ind_apices[i]}};
    for (int s= 0; s < 3; s++)
    {
      triangle  cur_facet = triangle_from_index_cube(facets[s], list->dim);
      tri_index idx;
      indices_unique_cube(facets[s][0], facets[s][1], facets[s][2], idx); //Convert triangle to it's data_structure index
      /*
       * Insert all of these sides into the tri_list. Note that this may cause race conditions!
       * We avoid these by looking up the first vertex of this triangle, and aquiring this lock.
       */  
      omp_set_lock(&locks[idx[0]][idx[1]]);//Thread-safe now! As we edit the data check_list.p_arr[idx[0]][idx[1]]
      if (tri_list_insert(list, &cur_facet, TRI_LIST_NO_RESIZE))
        facets_add_cnt++;
      omp_unset_lock(&locks[idx[0]][idx[1]]);
    }
  }
}
/*
 * This function removes all the triangles from the data list that intersect with the given tetra.
 * All the possible new triangles that might be non-conform are added to the check_list.
 */
size_t filter_intersection_data_list_tet(data_list * data, tri_list * check_list, ptetra tet, omp_lock_t ** locks) {
  int dim = data_list_dim(data);
  size_t result = 0;

  cube_points cube = gen_cube_points(dim);
  static facet_acute_data parameters;
  #pragma omp threadprivate(parameters)
  #pragma omp parallel 
  {
    parameters.cube = &cube;
    parameters.boundary_func = &triangle_boundary_cube;
    parameters.data = data;
    parameters.store_acute_ind = 1;
    parameters.acute_ind  = malloc(sizeof(vert_index) * cube.len);
  }

  if (tri_list_count(check_list) > 0)
    printf("The check list is not empty.. Oke?");

  tri_index cur_idx;
  triangle cur_tri;
  size_t i,j;
  int k;
  #pragma omp parallel for  schedule(dynamic,dim) shared(locks) private(cur_tri, cur_idx, i,j,k)
  for (i = 0; i < cube.len; i++) 
    for (j = 0; j < cube.len - i; j++)
    {
      //if mem_list->t_arr[i][j] does not exist
      if (data->mode == DATA_MEM_LIST_CUBE && !data->mem_list.t_arr[i][j])
        continue;

      for (k = data_list_dim_size(data, 2, i, j) - 1; k >= 0; k--) {
        cur_idx[0] = i;
        cur_idx[1] = j;
        if (data->mode == DATA_MEM_LIST_CUBE)
        {
          cur_idx[2] = k;
          if (!GMI(data->mem_list.t_arr, cur_idx))
            continue; 
        } else 
          cur_idx[2] = data->list.t_arr[i][j].p_arr[k];

        cur_idx[1] += i;
        cur_idx[2] += cur_idx[1];

        cur_tri = triangle_from_index_cube(cur_idx, dim);
        if (!tri_tet_disjoint(&cur_tri, tet)){ //If not disjoint with new tetrahedron, delete
          facet_conform(&cur_tri, &parameters);
          facets_tetra_list(check_list, cur_idx, parameters.acute_ind, parameters.acute_ind_len, locks);
          //Add all the sides of conform tetrahedrons with cur_tri as base to the possible non-conform list.
          //Cur_tri is not conform, remove from the data structure.
          if (data->mode == DATA_MEM_LIST_CUBE) {
            mem_list_cube_clear(&data->mem_list, &cur_tri);
          } else {
            tri_list_remove(&data->list, &cur_tri, TRI_LIST_NO_RESIZE);
          }

          #pragma omp atomic
          result++;
        }
      }
    }
  #pragma omp parallel 
  {
    free(parameters.acute_ind);
  }
  free(cube.points);
  return result;
}

void facets_conform_dynamic_remove(data_list * data, ptriangulation triang, int iterations, tri_list * check_list, tri_list * check_list_new, omp_lock_t ** locks) {
  int dim = data_list_dim(data);
  //tri_mem_list * list = &data->mem_list;

  cube_points cube = gen_cube_points(dim);
  //Initalize the parameter. Every thread should have it's own copy
  static facet_acute_data parameters;
#pragma omp threadprivate(parameters)
#pragma omp parallel 
  {
    parameters.cube = &cube;
    parameters.boundary_func = &triangle_boundary_cube;
    parameters.data = data;
    parameters.store_acute_ind = 0;
    parameters.acute_ind  = malloc(sizeof(vert_index) * cube.len);
  }


  int iter = 0;
  double time_start, time_check;
  size_t count = 0;
  int triang_consistent = 1;
  while (tri_list_count(check_list) && triang_consistent && (iter != iterations)) //While we have triangles to be removed)
  {
    time_start = omp_get_wtime();
    triangle cur_tri;
    tri_index cur_idx;
    int l,k;
    size_t i,j;
    size_t facets_add_total = 0;


#pragma omp parallel shared(locks) private(cur_tri, cur_idx, i,j,k,l)
    {
      facets_add_cnt = 0;
      if (omp_get_thread_num() == 0) {
        size_t new_count = data_list_count(data);
        if (count)
          printf("Removed %zu triangles\n", count - new_count);
        printf("\n\nLoop %d of conform dynamic\n", iter++);
        printf("Size of entire list    %zu\n", new_count);
        printf("Size of check list     %zu\n", tri_list_count(check_list));
        tri_list_validate(check_list);
        count = new_count;
      }
      /*
       * Loop over all the triangles in the check list. Check if they are not conform
       * if so, add all the possible new non-conform edges to the tmp_check_list.
       */
      #pragma omp for  schedule(dynamic,dim) 
      for (i = 0; i < cube.len; i++) {
        if (!triang_consistent) 
          continue; //We want break, but that is not possible with openMP

        for (j = i; j < cube.len; j++) {
          if (!triang_consistent)
            break;
          for (l = check_list->t_arr[i][j- i].len - 1; l >= 0; l--) {  //Loop over all triangles (i,j,*)
            k = check_list->t_arr[i][j - i].p_arr[l] + j;
            cur_idx[0] = i;
            cur_idx[1] = j;
            cur_idx[2] = k;
            cur_tri = triangle_from_index_cube(cur_idx, dim);
            //Cur_tri now holds the triangle we should check
            if ((data->mode == DATA_TRI_LIST && !tri_list_contains(&data->list,&cur_tri)) ||
                (data->mode == DATA_MEM_LIST_CUBE && !mem_list_cube_contains(&data->mem_list,&cur_tri)))
              continue; //This triangle was already removed.. Skip :-)

            if (!facet_conform(&cur_tri, &parameters)) { //This triangle is not conform, delete!
              parameters.store_acute_ind = 1;
              facet_conform(&cur_tri, &parameters);
              parameters.store_acute_ind = 0;
              //Add all the sides of conform tetrahedrons with cur_tri as base to the possible non-conform list.
              facets_tetra_list(check_list_new, cur_idx, parameters.acute_ind, parameters.acute_ind_len, locks);
              //Cur_tri is not conform, remove from the data structure.
              if (data->mode == DATA_MEM_LIST_CUBE)
                mem_list_cube_clear(&data->mem_list, &cur_tri);
              else
                tri_list_remove(&data->list, &cur_tri, TRI_LIST_NO_RESIZE);
            }
          }
        }

        if (omp_get_thread_num() == 0)
          triang_consistent = triangulation_consistent(triang, &parameters);
      }

      #pragma omp atomic
      facets_add_total += facets_add_cnt;
    }

    if (triang_consistent) {
      printf("Amount of triangles in new_check_list: +/-   %zu\n", facets_add_total);
      printf("Amount of triangles in new_check_list: exact %zu\n", tri_list_count(check_list_new));
      //Checked all the triangles from check_list. Empty it and swap the lists.
      tri_list_empty(check_list);

      tri_list tmp = *check_list;
      *check_list = *check_list_new;
      *check_list_new = tmp;
    } else
      printf("Triangulation not consistent anymore\n");

    time_check = omp_get_wtime();
    printf("\nTook %f seconds to construct new check list\n",time_check - time_start);
  }

  free(cube.points);
  #pragma omp parallel
  {
    free(parameters.acute_ind);
  }
}

ptriangle triangulation_start_facet(data_list  * data) {
  int dim = data_list_dim(data);
  /* Try to find a start facet on the z = 0 plane */
  ptriangle  start_facet = calloc(sizeof(triangle), 1);
  int found_start = 0;
  while (!found_start) { 
    /*
     * Random triangle with vertices:
     * (0,0,0), (?,0,0), (?,?,0)
     */
    start_facet->vertices[1][0] = rand() %dim;
    start_facet->vertices[2][0] = rand() %dim;
    start_facet->vertices[2][1] = rand() %dim;

    if (!triangle_acute(start_facet))
      continue;
    if (data->mode == DATA_MEM_LIST_CUBE)
      found_start = mem_list_cube_contains(&data->mem_list,start_facet);
    else
      found_start = tri_list_contains(&data->list, start_facet);
  }
  return start_facet;
}

triangulation triangulate_cube_random(data_list * data) {
  int dim = data_list_dim(data);
  cube_points cube = gen_cube_points(dim); 
  triangulation result;
  omp_lock_t result_lock;
  omp_init_lock(&result_lock);  //If we found a triangulation, use this lock!

  facet_acute_data parameters; //Parameters for conform_
  triangulation tmp_triang; //Triangulation we are expanding in current thread
  ptetra tet_list; //List of tetrahedrons, used in the parallel section
  unsigned short tet_list_len; //Holds the length of this list
  int triangulation_found = 0; //Stop if one of the threads has found a triangulation!
  int rand_bound, i;
  unsigned short tet_max, tet_min, tet_rand, tet_add;

  size_t max_volume;
  //Start the parallel loop!
#pragma omp parallel default(none) \
  private(parameters, tmp_triang, tet_list, tet_list_len, rand_bound, i,max_volume,tet_max, tet_min, tet_rand, tet_add) \
  shared(result, result_lock, cube,data,dim, triangulation_found)
  {
    //Initalization for each thread
    parameters.cube = &cube;
    parameters.boundary_func = &triangle_boundary_cube;
    parameters.data = data;
    parameters.store_acute_ind = 1;
    parameters.acute_ind  = malloc(sizeof(vert_index) * cube.len);

    tet_list = malloc(sizeof(tetra) * cube.len);
    max_volume = 0;

    while (!triangulation_found) { //Not found a triangulation
      //Initalize the triangulation variables
      tmp_triang = triangulation_init(dim);
      tet_list_len = 0;
      //Start triangle (0,0,0), (rand,0,0), (rand,rand,0)
      tmp_triang.bound_len = 1;
      tmp_triang.bound_tri = triangulation_start_facet(data);
      //printf("Thread %d with iteration %zu starts with:\n", omp_get_thread_num(), ++iterations);
      //print_triangle(tmp_triang.bound_tri);

      //While we have triangles on the boundary
      while (tmp_triang.bound_len > 0) {
        /*
         * We are going to add a tetrahedron on the boundary triangle.
         * To do so, we select a random triangle on the boundary. Then we generate all the
         * acute tetrahedra (above and below) with facets in our possible list.
         * From this list we remove all the tetrahedrons that intersect with our current triangulation.
         * Then we add a random tetrahedron to our triangulation and repeat.
         */
        rand_bound = rand() % tmp_triang.bound_len;
        //
        //Calculate the conform tetrahedrons above and below
        if (!facet_conform(tmp_triang.bound_tri + rand_bound, &parameters))
          break; //Triangle on the boundary that does not have a conform facet

        tet_list_len = parameters.acute_ind_len;
        //Form explicit list of the tetrahedrons
        for (i = 0; i < tet_list_len; i++) 
        {
          copyArr3(tet_list[i].vertices[0], tmp_triang.bound_tri[rand_bound].vertices[0]);
          copyArr3(tet_list[i].vertices[1], tmp_triang.bound_tri[rand_bound].vertices[1]);
          copyArr3(tet_list[i].vertices[2], tmp_triang.bound_tri[rand_bound].vertices[2]);
          copyArr3(tet_list[i].vertices[3], cube.points[parameters.acute_ind[i]]);
        }

        //Remove all the tetrahedrons that intersect with current triangulation.
        filter_tet_list_disjoint_triangulation(tet_list, &tet_list_len, &tmp_triang);

        if (tet_list_len == 0) 
          break; //We can not find a conform tetrahedron for this boundary.. Restart

        //Select a ttetrahedron from the tet_list to add to the triangulation.. Different approaches.
        //Combinations between: random tetra, smallest volume, maximum volume. Indices stored in tet_max, tet_min and tet_rand

        tet_list_min_max_volume(tet_list, tet_list_len, &tet_max, &tet_min);
        tet_rand = rand() % tet_list_len; 

        switch (omp_get_thread_num() % 6) {
          case 0: //Choose tet with max volume
            tet_add = tet_max; break;
          case 1: //Choose tet with min volume
            tet_add = tet_min; break;
          case 2: //Choose random tet
            tet_add = tet_rand; break;
          case 3: //Choose either max or min (random)
            tet_add = (rand() % 2)? tet_min : tet_max; break;
          case 4: //Choose either max or rand
            tet_add = (rand() % 5)? tet_max : tet_rand; break;
          case 5: //Either min or rand
            tet_add = (rand() % 5)? tet_min : tet_rand; break;
          default:
            tet_add = 0;
        }

        /*
         * Add the above tetra to the triangulation.
         * This removes all the boundary triangles that are covered by this tetrahedron
         */
        add_tet_triangulation(tet_list + tet_add,&tmp_triang);
      }
      if (triangulation_volume(&tmp_triang) > max_volume) {
        max_volume = triangulation_volume(&tmp_triang);
        printf("Record for thread %d using method %d amount: %zu\n", omp_get_thread_num(), omp_get_thread_num() % 6, max_volume);
        triangulation_print(&tmp_triang);
      }
      if (tmp_triang.bound_len == 0)
      {
        printf("FOUND A TRIANGULATION!!!\n");
        triangulation_print(&tmp_triang);
        omp_set_lock(&result_lock);
        result = tmp_triang;
        triangulation_found = 1;
        omp_unset_lock(&result_lock);
      } else 
        triangulation_free(&tmp_triang);
    }
    free(parameters.acute_ind);
    free(tet_list);
  }
  free(cube.points);
  omp_destroy_lock(&result_lock);
  return result;
}



triangulation triangulate_cube(data_list * data,  char * tmp_triang_file, char * tmp_data_file) {
  triangulation result = triangulation_init(data_list_dim(data));

  cube_points cube = gen_cube_points(result.dim);
  facet_acute_data parameters;
  parameters.cube = &cube;
  parameters.boundary_func = &triangle_boundary_cube;
  parameters.data = data;
  parameters.store_acute_ind = 1;
  parameters.acute_ind = malloc(sizeof(unsigned short) * cube.len);

  //This list holds all conform tetrahedrons for a given triangle, max size = cube.len
  ptetra tet_list = malloc(sizeof(tetra) * cube.len);
  unsigned short tet_list_len = 0;

  //Lists needed for the dynamic_remove loop
  tri_list check_list, check_list_new;

  check_list     = tri_list_init(result.dim, MEM_LIST_FALSE);
  check_list_new = tri_list_init(result.dim, MEM_LIST_FALSE);

  //Start triangle (0,0,0), (rand,0,0), (rand,rand,0)
  result.bound_len = 1;
  result.bound_tri = triangulation_start_facet(data);
  printf("Starting triangulation with facet:\n");
  print_triangle(result.bound_tri);
  /*
   * During this method we are going to operate data that is not thread-safe.
   * To avoid race conditions we need an array of locks. We use a lock for the
   * first two points of a triangle (so need 2d array of locks).
   */
  omp_lock_t ** locks = malloc(sizeof(omp_lock_t *) * cube.len);
  //Initalize the locks
  for (size_t i = 0; i < cube.len; i++){
    locks[i] = malloc(sizeof(omp_lock_t) * (cube.len - i));
    for (size_t j = 0; j < cube.len - i; j++)
      omp_init_lock(&locks[i][j]);
  }
  //While we have triangles on the boundary..
  while (result.bound_len > 0) {
    tri_list_empty(&check_list);
    tri_list_empty(&check_list_new);
    /*
     * We are going to add a tetrahedron on the boundary triangle.
     * To do so, we select a random triangle on the boundary. Then we generate all the
     * acute tetrahedra (above and below) with facets in our possible list.
     * From this list we remove all the tetrahedrons that intersect with our current triangulation.
     * Then we add a random tetrahedron to our triangulation, update the conform list and repeat.
     */
    int rand_bound = rand() % result.bound_len;
    printf("\n\nTotal amount of triangles left:%zu\nExpanding triangulation at boundary triangle: \n", data_list_count(data));
    print_triangle(result.bound_tri + rand_bound);

    //Calculate the conform tetrahedrons above and below
    if (!facet_conform(&result.bound_tri[rand_bound], &parameters))
    {
      printf("We have a triangle on the boundary that is not conform anymore.\n");
      printf("Whatthefuck? Breaking!\n");
      break;
    }

    tet_list_len = parameters.acute_ind_len;
    printf("Total amount of conform tetrahedrons found for this boundary: %hu\n", tet_list_len);
    //Form explicit list of the tetrahedrons
    for (unsigned short i = 0; i < tet_list_len; i++) 
    {
      copyArr3(tet_list[i].vertices[0], result.bound_tri[rand_bound].vertices[0]);
      copyArr3(tet_list[i].vertices[1], result.bound_tri[rand_bound].vertices[1]);
      copyArr3(tet_list[i].vertices[2], result.bound_tri[rand_bound].vertices[2]);
      copyArr3(tet_list[i].vertices[3], cube.points[parameters.acute_ind[i]]);
    }

    //Remove all the tetrahedrons that intersect with current triangulation.
    filter_tet_list_disjoint_triangulation(tet_list, &tet_list_len, &result);

    printf("Amount of tetrahedrons left after filtering: %hu\n\n",tet_list_len);
    if (tet_list_len == 0) {
      printf("Waarom is deze lijst nu al fucking leeggefilterd?\n");
      printf("Dead end, helaas pindakaas. Got to %zu\n", result.tetra_len);
      break;
    }

    //Select random tetrahedron disjoint with the current triangulation
    int rand_tet = rand() % tet_list_len;
    /*
     * Add the above tetra to the triangulation.
     * This removes all the boundary triangles that are covered by this tetrahedron
     */
    printf("Adding the following tetra to the triangulation\n");
    print_tetra(tet_list + rand_tet);
    printf("\n\n");
    add_tet_triangulation(tet_list + rand_tet, &result);
    triangulation_print(&result);

    if (!result.bound_len) //If we have no boundaries left, we must be done!!
    {
      printf("No more boundaries left.. WE FINNISHED!??\n");
      break;
    }
    //Consistency check
    if (!triangulation_consistent(&result, &parameters))
    {
      printf("Triangulation not consistent after adding the tetrahedron. Breaking.\n");
      break;
    }
    /*
     * Calculate a list of all the triangles we are going to remove
     */
    size_t removed = filter_intersection_data_list_tet(data,  &check_list, tet_list + rand_tet, locks);
    printf("Removed %zu triangles that are not disjoint with the new tetrahedron\n", removed);
    printf("The check_list has size %zu\n", tri_list_count(&check_list));

    if (!triangulation_consistent(&result, &parameters)) {
      printf("After filtering the memory list we have a non consistent triangulation. Break\n");
      break;
    }
    //Do two iterations
    facets_conform_dynamic_remove(data, &result, 1, &check_list, &check_list_new, locks);

    if (!triangulation_consistent(&result, &parameters)) {
      printf("Triangulation not consistent anymore after conforming the data set.. Breaking\n");
      break;
    }

    /*mem_list_cube_compress(&data->mem_list);


      if (tmp_triang_file && tmp_data_file) {
      triangulation_to_file(&result, tmp_triang_file);
      data_list_to_file(data, tmp_data_file, MEM_LIST_SAVE_CLEAN);
      }
      */
  }
  for (size_t i = 0; i < cube.len; i++){
    for (size_t j = 0; j < cube.len - i; j++)
      omp_destroy_lock(&locks[i][j]);
    free(locks[i]);
  }

  free(locks);
  free(cube.points);
  free(parameters.acute_ind);
  free(tet_list);
  tri_list_free(&check_list);
  tri_list_free(&check_list_new);
  printf("Triangulation has length of %zu\n", result.tetra_len);
  return result;
}
/*
   void filter_tri_list_remove_list(tri_list * list, tri_list * remove_list) {
   size_t dim_size = tri_list_dim_size(list->dim);
   triangle cur_tri;
   size_t i,j,k;
   int l;

#pragma omp parallel for schedule(dynamic,list->dim) private(j,k,i,l,cur_tri) 
for (i = 0; i < dim_size; i++) {
vertex_from_index_cube(i,remove_list->dim, cur_tri.vertices[0]);
for (j = i; j < dim_size; j++) {
vertex_from_index_cube(j, remove_list->dim, cur_tri.vertices[1]);
for (l = remove_list->t_arr[i][j-i].len - 1; l >= 0; l--) {  //Loop over all triangles (i,j,*)
k = remove_list->t_arr[i][j-i].p_arr[l] +  j;
vertex_from_index_cube(k, remove_list->dim, cur_tri.vertices[2]);
//Cur_tri is the l-th triangle with points(i,j,*)
//Remove cur_tri from list
tri_list_remove(list, &cur_tri, TRI_LIST_NO_RESIZE);
}
}
}
}
*/
/*
 *
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
