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
 * All normals point inward or outward... PROOF
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
  if (dotArr3(P[0],P[1]) < 0 || dotArr3(P[0], P[2]) < 0 || dotArr3(P[1],P[2]) < 0 )
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

//For history sake!
int tetra_acute_old(ptetra tet) {
  //First check if it's facets are acute.
  //Old code
  /*
  if (!mat3_triangle_acute(tet->vertices) || //facet 0,1,2
      !mat3_triangle_acute(&tet->vertices[1]) || //facet 1,2,3
      !arr3_triangle_acute(tet->vertices[0],tet->vertices[1],tet->vertices[3]) || //facet 0,1,3
      !arr3_triangle_acute(tet->vertices[0],tet->vertices[2],tet->vertices[3])) //facet 0,2,3
    return 0;
  mat3 P;
  subArr3(tet->vertices[1], tet->vertices[0], P[0]);
  subArr3(tet->vertices[2], tet->vertices[0], P[1]);
  subArr3(tet->vertices[3], tet->vertices[0], P[2]);
  */
  //Optimized "hacky-code". Does the same as above, but then entirely inlined. As this function is called the most, this 
  //means some performance gain!
  
  arr3 P[6]; //All the edges
  
  #define triangle_side(A,B,P) subArr3(tet->vertices[A],tet->vertices[B],P)
  
  
   // Calculate the sides of our tetrahedron check if they are acute
   
   
  //Facet 0,1,2
  triangle_side(1,0,P[0]); //P[0] = v1 - v0
  triangle_side(2,0,P[1]); //P[1] = v2 - v0
  triangle_side(2,1,P[3]); //P[3] = v2 - v1 NOTE P[3], P[2] is saved for something else
  //Convention third argument must be the side vector from v1 to v2
  if (!triangle_sides_acute(P[0],P[1],P[3])) 
    return 0;
    
  //Facet 0,1,3  
  triangle_side(3,0,P[2]); //P[2] = v3 - v0
  triangle_side(3,1,P[4]); //P[4] = v3 - v1
  if (!triangle_sides_acute(P[0], P[2], P[4]))
    return 0;
  
  //Facet 0,2,3
  triangle_side(3,2,P[5]); //P[5] = v3 - v2
  if (!triangle_sides_acute(P[1], P[2], P[5]))
    return 0;
  //Facet 1,2,3
  if (!triangle_sides_acute(P[3], P[4], P[5]))
    return 0;   
  
  
  /*
   * Get first three sides. P already holds V1-V0 and V2-V0. Outward pointing normal
   */
  //subArr3(tet->vertices[1], tet->vertices[0], P[0]);
  //subArr3(tet->vertices[2], tet->vertices[0], P[1]);
  //subArr3(tet->vertices[3], tet->vertices[0], P[2]);
  mat3 Q_cross;
  crossArr3(P[1],P[2], Q_cross[0]);
  crossArr3(P[2],P[0], Q_cross[1]);
  crossArr3(P[0],P[1], Q_cross[2]);
  
  //Normals on three facets, point outward
  //q1 loodrecht op p1, p2
  //als q1 ook loodrecht op p0, dan liggen drie vectors in hetzelfde vlak, degenerate..
  if (dotArr3(P[0], Q_cross[0]) == 0)
    return 0;
  //Calculate inner products <q_i, q_j>
  mat3 Q_cross_trans;
  transMat3(Q_cross, Q_cross_trans);
  //Q_cross has normals in rows, so transpose the right one to get inner product matrix
  mat3 Q_prod;
  dotMat3(Q_cross,Q_cross_trans, Q_prod);
  
  //Check of eerste drie facetten niet-scherpe hoek maken
    if (Q_prod[0][1] >= 0 || Q_prod[0][2] >= 0 || Q_prod[1][2] >= 0)
      return 0;
    
  //Controleer of de andere drie facetten een niet-scherpe hoek maken
  vec3 Q_last = mat3Vector(Q_prod, (vec3) {{-1,-1,-1}});
  return (Q_last.pt[0] < 0 && Q_last.pt[1] < 0 && Q_last.pt[2] < 0);
}

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

/*
 * Returns whether this tetrahedron only has facets in the mem_list
 */
 
int tetra_facets_list(ptriangle triang, arr3 new_vertex, tri_mem_list * acute_list) {
  //triang itself is included in the list.. Otherwise dumb call!
  tri_index indices[3];
  
  if (acute_list->fund) //We need to symmetry the triangle to a fundamental-area representation first
    return (mem_list_get_sym_fund(acute_list, triang->vertices[0], triang->vertices[1], new_vertex) &&
            mem_list_get_sym_fund(acute_list, triang->vertices[0], triang->vertices[2], new_vertex) &&
            mem_list_get_sym_fund(acute_list, triang->vertices[1], triang->vertices[2], new_vertex)); 
  
  else {
    vertices_to_index(triang->vertices[0], triang->vertices[1], new_vertex, acute_list->dim_mult, indices[0]);
    vertices_to_index(triang->vertices[0], triang->vertices[2], new_vertex, acute_list->dim_mult, indices[1]);
    vertices_to_index(triang->vertices[1], triang->vertices[2], new_vertex, acute_list->dim_mult, indices[2]);
    for (int i = 0; i < 3; i++) {
      if (!EMI(acute_list->t_arr, indices[i]) ||  //Point does not exist or is not acute
          !GMI(acute_list->t_arr, indices[i]))
        return 0;
     } 
    return 1;
  }
}



/*
 * Returns whether this triangle has an acute tetrahedron above and below. More information
 * is stored in result. Convention: normal points to above (inwards for triangles on the boundary).
 * 
 * If acute_list is set, only tetrahedrons with facets in this list may be included.
 */
int triangle_tetra_acute(ptriangle triang, cube_points * cube, ptriang_tetra_result res  , tri_mem_list *acute_list ) {
  //Check if triangle is acute..
  if (!mat3_triangle_acute(triang->vertices)) //Can replace with code below
    return 0;
    
 /* 
  * Initalize variables
  */
  res->boundary_above = 0;
  res->boundary_below = 0;
  res->acute_above    = 0;
  res->acute_below    = 0;
  arr3 * cube_pts = cube->points;
    
  mat3 P;
  //Store sides in matrix P
  //P[0] = V1 - V0
  //P[1] = V2 - V0
  //P[2] = V2 - V1
  triangle_sides(triang->vertices[0], triang->vertices[1], triang->vertices[2],P);
  
  /*
   * Calculate the normals for this triangle
   * - The normal on the triangle plane, to determine on which side a point is
   * - The three normals on the "side"-planes, used to check if a point lies in the toblerone
   *
   * Now calculate the plane equations normal . p = d
   */
  //Convention P[1] cross P[0]
  arr3 tri_normal;
  crossArr3(P[1], P[0], tri_normal);
  int  tri_d = dotArr3(tri_normal, triang->vertices[0]); //For the plane equation of the triangle
  //triangle_normal(triang, normal);
  
  
  mat3 side_normals;
  arr3 side_d; //The constant expression for the plane equations of the sides
  //Convention, need to explain why it works?? Third point must be on the other side!!
  
  crossArr3(tri_normal, P[0], side_normals[0]);
  crossArr3(P[1], tri_normal, side_normals[1]);
  crossArr3(tri_normal, P[2], side_normals[2]);
  
  side_d[0] = dotArr3(side_normals[0], triang->vertices[0]); 
  side_d[1] = dotArr3(side_normals[1], triang->vertices[0]);
  side_d[2] = dotArr3(side_normals[2], triang->vertices[2]);
  
  int dotprod;
  //Determine whether triangle lies entirely on boundary plane.
  if (triangle_boundary(triang,cube->dim)) {
    arr3 testpt = {1,1,1}; //Point inside the cube
    /*
     * Calculate angle between normal and vector from vertex to point inside cube (testpt)
     * If both point the same way, the angle < 90 degrees, and the dotproduct must be positive.
     */
    dotprod = dotArr3(testpt, tri_normal);
    if (dotprod > tri_d) {
    /*
     * Normal points inward. So the boundary lies "below" the triangle. And therefore we don't need to check below
     */
      res->boundary_below = 1;
      res->acute_below    = 1;
    } else {
    /*
     * Normal points outside the cube. Boundary lies "above". We only need to check "below" the triangle"
     */
      res->boundary_above  = 1;
      res->acute_above     = 1;
    } 
  }
  
  for (size_t i = 0; i < cube->len; i++){
    dotprod = dotArr3(cube_pts[i], tri_normal); //Calculate angle between normal and vector from plane to cube_pt
    
    
    //Check if current point lies on side of the triangle which isn't sharp yet and check if point lies in the
    //prism of possible sharp tetrahedron
    
    if (((dotprod > tri_d && !res->acute_above) ||  //Dotprod > tri_d implies that the point lies "above" the triangle (same side as normal) 
         (dotprod < tri_d && !res->acute_below))     //Dotprod < tri_d implies lies below
          && //Check if point lies in the prism
          dotArr3(cube_pts[i],side_normals[0]) < side_d[0] &&
          dotArr3(cube_pts[i],side_normals[1]) < side_d[1] &&
          dotArr3(cube_pts[i],side_normals[2]) < side_d[2]) 
      {
       if (acute_list && !tetra_facets_list(triang, cube_pts[i], acute_list))
            continue;
      tetra test_tetra;
      memcpy(test_tetra.vertices, &cube_pts[i], sizeof(arr3));
      memcpy(test_tetra.vertices + 1, triang->vertices, 3 * sizeof(arr3)); 
      //test_tetra now holds the four vertices.
      if (tetra_acute(&test_tetra)) {//Yes this one is acute!
        if (dotprod > tri_d) {
          res->acute_above = 1;
        } else {
          res->acute_below  = 1;
        }
        if (res->acute_above && res->acute_below) {
          return 1;
        }
      }
    }
  }
  return 0;  
}



/*
 * Dindex * comb_list = combinations_list(cube.len, 2,&len); //All possible pairs
    cube_points fund = gen_fund_points(dim); 
    index_list = malloc(sizeof(arr3) * len * fund.len);
    c = 0;
    for (size_t i = 0; i < len; i++) {
      triangle cur_triangle;
      Dindex v1,v2;
      v1 = comb_list[i * 2];
      v2 = comb_list[i * 2 + 1];
      memcpy(cur_triangle.vertices, cube.points[v1], sizeof(arr3));
      memcpy(cur_triangle.vertices + 1, cube.points[v2], sizeof(arr3));
      //First two points of cur_triangle now hold the combination index i
      for (int j = 0; j < fund.len; j++)
      {
        memcpy(cur_triangle.vertices + 2, fund.points[j], sizeof(arr3));
        arr3 index;
        triangle_to_index(&cur_triangle, acute_list->dim, index);
        if (EMI(acute_list->t_arr,index) && GMI(acute_list->t_arr,index)) {//Triangle is acute
          memcpy(index_list + c, index, sizeof(arr3));
          c++;
        }
      }
    }
    free(comb_list);

*/


/*
 * Filter list on face2face property. 
 * If sym is true, it only checks for triangles in the fundamental domain (and then filters the 48 symmetries).
 * If sym is false, it checks all triangles that are set in the acute_list
 * BROKEN AS FUCK
 */
void mem_list_face2face_old(tri_mem_list * acute_list, int sym){
  cube_points cube = gen_cube_points(acute_list->dim); 
  
  int *  fund_index = NULL;
  size_t fund_len   = 0;
  //Calculate the indices of the fundamental domain
  if (sym) {
    cube_points fund = gen_fund_points(acute_list->dim[0]);
    fund_len = fund.len;
    fund_index = malloc(sizeof(int) * fund_len);
    for (size_t i = 0; i < fund_len; i++) 
      fund_index[i] = vertex_to_index(fund.points[i], acute_list->dim_mult);
    free(fund.points);
  }
  
  int changed = 1;
  triang_tetra_result tri_result;
  while (changed) {
    changed = 0;
    printf("Start filter loop, amount of acute facets: %zu\n", mem_list_count(acute_list));
    //Loop over all indices
    int vertex_fund = 0; //Indicates whether we have a vertex in the fundamental domain
    for (int i = 0; i < acute_list->dim_size[0]; i++) {
      if (!acute_list->t_arr[i])
        continue;
      if (sym)
        vertex_fund = index_in_array(i, fund_index, fund_len);
      for (int j = 0; j < acute_list->dim_size[1]; j++) { 
        if (acute_list->t_arr[i][j]) { //Third dimension exists!
          if (sym && !vertex_fund)
            vertex_fund = index_in_array(j, fund_index, fund_len);
          if (sym && !vertex_fund) // Only loop over fundamental domain!!
            for (size_t k = 0; k < fund_len; k++) {
              tri_index index = {i,j,fund_index[k]};
              if (!GMI(acute_list->t_arr, index))
                continue; //Triangle no longer in the acute list, no need to check
              triangle cur_tri = triangle_from_index(index, cube.dim);
              if (!triangle_tetra_acute(&cur_tri,&cube, &tri_result,acute_list)) { //This triangle is not acute anymore
                mem_list_clear_sym(acute_list, &cur_tri);
                changed = 1;
              }
            }
          else //Loop over entire domain
            for (int k = 0; k < acute_list->dim_size[2]*8; k++) {
              tri_index index = {i,j,k};
              if (!GMI(acute_list->t_arr, index))
                continue; //Triangle no longer in the acute list, no need to check
              triangle cur_tri = triangle_from_index(index, cube.dim);
              if (!triangle_tetra_acute(&cur_tri,&cube, &tri_result,acute_list)) { //This triangle is not acute anymore
                mem_list_clear_sym(acute_list, &cur_tri);
                changed = 1;
              }
            }
        }
      }
    }
  }
  if (sym)
    free(fund_index);
  free(cube.points);
}

void mem_list_face2face(tri_mem_list * acute_list){
  size_t len = 0;
  arr3 * cube_pts;
  
  
  cube_points fund_points = gen_fund_points(acute_list->dim[0]);
  arr3 *fund_pts = fund_points.points;
  cube_points cube = gen_cube_points(acute_list->dim);
  len = cube.len;
  cube_pts = cube.points;
  Dindex	*arr;
  int changed = 1;
  while (changed) {
    changed = 0;
    printf("Starting face2face loop with %zu acute triangles\n", mem_list_count(acute_list));
    //Check every point in the fundamental domain against all other points 
    arr = revdoor_init(len, 2);  
    if (arr == 0) {
      puts("Error involving revdoor_init");
      exit(1);
    }
    triang_tetra_result tri_result;
    //Loop through all combinations
    do {
      size_t j = arr[0], k = arr[1];
      tri_index indices;
      for (size_t i = 0; i < fund_points.len; i++){
        //Check if the current triangle is in the memory_list "acute"
        if (acute_list->fund) {
          //if (!mem_list_get_sym_fund(acute_list,fund_pts[i],cube_pts[j],cube_pts[k])) //Checks all symmetries
            //continue; 
          vertices_to_index_fund(fund_pts[i],cube_pts[j],cube_pts[k], acute_list, indices); 
        } else
          vertices_to_index(cube_pts[j],cube_pts[k], fund_pts[i], acute_list->dim_mult, indices);
        if (!EMI(acute_list->t_arr, indices) || !GMI(acute_list->t_arr,indices))
          continue;
          
        triangle cur_tri = (triangle) {{{fund_pts[i][0],fund_pts[i][1],fund_pts[i][2]},
                                  {cube_pts[j][0],cube_pts[j][1],cube_pts[j][2]},
                                  {cube_pts[k][0],cube_pts[k][1],cube_pts[k][2]}}};
        if (!triangle_tetra_acute(&cur_tri,&cube, &tri_result,acute_list)) { //remove from list
          changed = 1;
          if (acute_list->fund)
            mem_list_clear_sym_fund(acute_list, &cur_tri);
          else
            mem_list_clear_sym(acute_list, &cur_tri);

        }
      }
    } while (revdoor_next(arr));
    revdoor_free(arr);
  } 
  free(cube_pts);
  free(fund_pts);
}


tri_mem_list acute_triangles_tetra(arr3 dim) {
  if (dim[0] != dim[1] || dim[1] != dim[2]) {
    printf("Cannot execute acute_triangles_tetra with different dimensions\n");
    exit(0);
  }
  /*
  tri_mem_list result = {NULL, 0, {dim[0],dim[1],dim[2]}}; 
  result.index_list = malloc(1000 * sizeof(tri_index));
  result.len = 1000;
  size_t count = 0;
  */  
  tri_mem_list result = mem_list_init_fund(dim[0],MEM_LIST_FALSE);
  
  cube_points fund = gen_fund_points(dim[0]);
  cube_points cube = gen_cube_points(dim);
  
  Dindex	*arr;
  //Check every point in the fundamental domain against all other points 

  arr = revdoor_init(cube.len, 2);  
  if (arr == 0) {
    puts("Error involving revdoor_init");
    exit(1);
  }
  
  triang_tetra_result tri_result;
  do{  
    //Loop through all combinations
    size_t j = arr[0], k = arr[1];  
    for (size_t i = 0; i < fund.len; i++){

      triangle cur_tri = (triangle) {{{fund.points[i][0],fund.points[i][1],fund.points[i][2]},
                                {cube.points[j][0],cube.points[j][1],cube.points[j][2]},
                                {cube.points[k][0],cube.points[k][1],cube.points[k][2]}}};
      if (triangle_tetra_acute(&cur_tri,&cube, &tri_result,NULL)) {
        /*
        if (count == result.len) {
          result.index_list = realloc(result.index_list, (count + count/100) * sizeof(tri_index));
          result.len = (count + count / 100);
        }          
        if (result.index_list == NULL) {
          puts("Error allocating memory!!");
          exit(1);
        }
        triangle_to_index(&cur_tri, dim_mult, result.index_list[count]);
        count++;
        */
        mem_list_set_sym_fund(&result, &cur_tri);
        //mem_list_set_sym(&result,&cur_tri);
      }
    } 
  } while (revdoor_next(arr));
  revdoor_free(arr);
  free(cube.points);
  free(fund.points);
  //result.len = count;
  return result;
}


triangle_list acute_triangles_tetra_old(arr3 dim) {
  triangle_list result = {NULL, 0, {dim[0],dim[1],dim[2]}};
  size_t count = 0;
  ptriangle t_arr = NULL;
  size_t len = 0;
  arr3 *cube_pts;
  
  cube_points cube = gen_cube_points(dim);
  len = cube.len;
  cube_pts = cube.points;
  
  
  Dindex	*arr;
   
  arr = revdoor_init(len, 3);  
  if (arr == 0) {
    puts("Error involving revdoor_init");
    exit(1);
    return result;
  }
  triang_tetra_result tri_result;
  
  //Loop through all combinations
  do {
    size_t j = arr[0], k = arr[1], i = arr[2];  
    triangle cur_tri = (triangle) {{{cube_pts[i][0],cube_pts[i][1],cube_pts[i][2]},
                              {cube_pts[j][0],cube_pts[j][1],cube_pts[j][2]},
                              {cube_pts[k][0],cube_pts[k][1],cube_pts[k][2]}}};
    if (triangle_tetra_acute(&cur_tri,&cube, &tri_result, NULL)) {
      
      count++;
      t_arr =  realloc(t_arr, count * sizeof(triangle));
      if (t_arr == NULL) {
        puts("Error allocating memory!!");
        exit(1);
      }
      t_arr[count - 1] = cur_tri;
    }
  } while (revdoor_next(arr));
  revdoor_free(arr);
  result.t_arr = t_arr;
  result.len = count;
  free(cube_pts);
  return result;
}

void print_tetra(ptetra tet) {
  printf("[[%d,%d,%d],\n",tet->vertices[0][0],tet->vertices[0][1],tet->vertices[0][2]);
  printf(" [%d,%d,%d],\n",tet->vertices[1][0],tet->vertices[1][1],tet->vertices[1][2]);
  printf(" [%d,%d,%d],\n",tet->vertices[2][0],tet->vertices[2][1],tet->vertices[2][2]);
  printf(" [%d,%d,%d]]\n",tet->vertices[3][0],tet->vertices[3][1],tet->vertices[3][2]);
}

/*
  cube_pts = [[x,y,z] for x in range(m+1) for y in range(n+1) for z in range(k+1)]
  combs = it.combinations(cube_pts,4)
  combs = np.array(list(combs))
  
  for comb in combs:
    tetra = Tetra(comb)
    acute = tetra.Acute()
    if acute:
      if (comb[:,1] == n).any():
        mat_edge = np.vstack((mat_edge,[comb]))
      mat = np.vstack((mat,[comb]))
  print mat_edge
  print mat.shape, mat_edge.shape
  return mat, mat_edge
  */





















int facet_tetra_list(ptriangle triang, arr3 new_vertex, tri_mem_list * acute_list) {
  return (mem_list_get_sym_fund(acute_list, triang->vertices[0], triang->vertices[1], new_vertex) &&
          mem_list_get_sym_fund(acute_list, triang->vertices[0], triang->vertices[2], new_vertex) &&
          mem_list_get_sym_fund(acute_list, triang->vertices[1], triang->vertices[2], new_vertex));
}

void tetra_add_array(tetra tetra_to_add, ptetra  * tetra_array, int * len) {
  (*len)++;
  *tetra_array = (ptetra) realloc(*tetra_array, (*len) * sizeof(tetra));
  if (*tetra_array == NULL) {
    puts("Error allocating memory!!");
    exit(1);
  }
  (*tetra_array)[*len - 1] = tetra_to_add;      
}

int facet_cube_acute(ptriangle triang, facet_acute_data * data, int mode) {
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
  
  data->boundary_triangle = triangle_boundary(triang,data->cube->dim); //Boundary plane only needs acute tetra on 1 side
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
        if ((mode == FACET_ACUTE_LIST) && !facet_tetra_list(triang, data->cube->points[i], data->acute_list))
          continue;   
        if (mode == FACET_ACUTE_TETRA) {//Explicitly create a list of tetrahedron
          if (dotprod > tri_d)
            tetra_add_array(test_tetra, &data->tetra_above, &data->tetra_above_len);
          else
            tetra_add_array(test_tetra, &data->tetra_below, &data->tetra_below_len);
        }
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
  #pragma omp parallel for schedule(guided) private(l,j,k,i,cur_tri) firstprivate(parameters)

  for (l = 0; l < comb_len; l++) {//Loop through all combinations
    j = comb[l*2];
    k = comb[l*2+1];
    for (i = 0; i < fund.len; i++){    
      cur_tri = (triangle) {{{fund.points[i][0],fund.points[i][1],fund.points[i][2]},
                             {cube.points[j][0],cube.points[j][1],cube.points[j][2]},
                             {cube.points[k][0],cube.points[k][1],cube.points[k][2]}}};
                      
      if (facet_cube_acute(&cur_tri,&parameters,FACET_ACUTE))
        mem_list_set_sym_fund(&result, &cur_tri);
    }
  }
  free(comb);
  free(cube.points);
  free(fund.points);
  return result;
}

//if save_file is set. The acute_list is saved to this file every hour
#define save_interval 30*60

void facets_face2face(tri_mem_list * acute_list, char * save_file){
  cube_points fund = gen_fund_points(acute_list->dim[0]);
  cube_points cube = gen_cube_points(acute_list->dim);
  size_t comb_len;
  Dindex * comb = combinations_list(cube.len, 2, &comb_len);
  int changed = 1;
  size_t l,i,j,k;
  tri_index indices;
  char tmp_file[100];
  if (save_file) 
    sprintf(tmp_file,"%s_tmp", save_file);
          
  unsigned short idx2, idx3;
  triangle cur_tri;
  facet_acute_data parameters;
  parameters.cube = &cube;
  parameters.acute_list = acute_list;
  double time_start =0 , time_end = 0, time_save = save_interval;
  while (changed) {
    changed = 0;
    printf("Face2face loop with %zu acute triangles from thread %d.\n"
          , mem_list_count(acute_list), omp_get_thread_num());
    time_start = omp_get_wtime();
    #pragma omp parallel for schedule(guided) private(l,j,k,i,cur_tri, idx2,idx3,indices)  firstprivate(parameters)
    for (l = 0; l < comb_len; l++) {//Loop through all combinations
      j = comb[l*2];
      k = comb[l*2+1];
      //OVERBODIGE STAP.. TOCH?
      idx2 = vertex_to_index(cube.points[j], acute_list->dim_mult);
      idx3 = vertex_to_index(cube.points[k], acute_list->dim_mult);
      for (i = 0; i < fund.len; i++){ //Against all fundamental points
        vertices_unique_fund(vertex_to_index(fund.points[i],acute_list->dim_mult), idx2, idx3, acute_list, indices); 
        if (!GMI(acute_list->t_arr,indices)) //Check if this index is still acute
          continue;
        cur_tri = (triangle) {{{fund.points[i][0],fund.points[i][1],fund.points[i][2]},
                             {cube.points[j][0],cube.points[j][1],cube.points[j][2]},
                             {cube.points[k][0],cube.points[k][1],cube.points[k][2]}}};
        if (!facet_cube_acute(&cur_tri,&parameters,FACET_ACUTE_LIST)) { //remove from list
          changed = 1;
          mem_list_clear_sym_fund(acute_list, &cur_tri);
        }
      }
      if (save_file && //Do we want to save the file?
          omp_get_thread_num() == 0 && //Only let master save to the file
        ((omp_get_wtime() - time_start) > time_save)) { //Time to save current progress
        
        printf("Saving tmp file with +/- %zu facets.\n", mem_list_count(acute_list));
        mem_list_to_file(acute_list, tmp_file, MEM_LIST_SAVE_CLEAN); //Save as tmp
        remove(save_file); //Remove the old progress file
        rename(tmp_file, save_file); //Rename tmp to new   
        time_save += save_interval;
      }
    } 
    time_end   = omp_get_wtime();
    printf("Loop took %f seconds.\n\n", time_end-time_start);
    
  }
  free(comb); 
  free(cube.points);
  free(fund.points);
}
