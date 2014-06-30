#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vector.h"
#include "triangle.h"
#include "combinations.h"

Triangle triangle_init(vec3 v1, vec3 v2, vec3 v3){
  Triangle result;
  result.vertices[0] = v1;
  result.vertices[1] = v2;
  result.vertices[2] = v3;
  /*
   * Holds the vectors of the three sides of the triangle
   * P[0] = AB, P[1] = AC, P[2] = BC
   */
  subArr3(v2.pt, v1.pt, result.P[0].pt);
  subArr3(v3.pt, v1.pt, result.P[1].pt);
  subArr3(v3.pt, v2.pt, result.P[2].pt);
  /*
  result.P[0] = subVectors(v2,v1);
  result.P[1] = subVectors(v3,v1);
  result.P[2] = subVectors(v3,v2);
  */
  return result;
}
//Vertices in rows of matrix
int mat3_triangle_acute(mat3 v){
  mat3 P;
  triangle_sides(v[0],v[1],v[2],P);
  return triangle_P_acute(P);
  /*
   * 
  subArr3(v[1],v[0], P[0]);
  subArr3(v[2],v[0], P[1]);
  subArr3(v[2],v[1], P[2]);
  return ((dotArr3(P[0],P[1]) > 0) && //Angle AB-AC
          (dotArr3(P[1],P[2]) > 0) && //Angle AC-BC=angle CA-CB
          ((-P[0][0] * P[2][0] - P[0][1]*P[2][1] - P[0][2]*P[2][2]) > 0)); 
        //  (dotArr3(negVector(P[0]),P[2]) > 0));  //Angle BA, BC
  */
}

int arr3_triangle_acute(arr3 v0, arr3 v1, arr3 v2) {
  mat3 P;
  triangle_sides(v0,v1,v2,P);
  return triangle_P_acute(P);
  /*
    subArr3(v1,v0, P[0]);
  subArr3(v2,v0, P[1]);
  subArr3(v2,v1, P[2]);
  return ((dotArr3(P[0],P[1]) > 0) && //Angle AB-AC
          (dotArr3(P[1],P[2]) > 0) && //Angle AC-BC=angle CA-CB
          ((-P[0][0] * P[2][0] - P[0][1]*P[2][1] - P[0][2]*P[2][2]) > 0)); 
        //  (dotArr3(negVector(P[0]),P[2]) > 0));  //Angle BA, BC  
  */
}
int triangle_acute(Triangle *triang) {
  
 return ((dotVector(triang->P[0],triang->P[1]) > 0) && //Angle AB-AC
        (dotVector(triang->P[1],triang->P[2]) > 0) && //Angle AC-BC=angle CA-CB
        (dotVector(negVector(triang->P[0]),triang->P[2]) > 0));  //Angle BA, BC
}

void triangle_normal(ptriangle triang, arr3 normal) {
  arr3 P[2];
  subArr3(triang->vertices[1],triang->vertices[0],P[0]);
  subArr3(triang->vertices[2],triang->vertices[0],P[1]);
  crossArr3(P[1], P[0], normal);
}

void print_triangle(ptriangle tet) {
  printf("[[%d,%d,%d],\n",tet->vertices[0][0],tet->vertices[0][1],tet->vertices[0][2]);
  printf(" [%d,%d,%d],\n",tet->vertices[1][0],tet->vertices[1][1],tet->vertices[1][2]);
  printf(" [%d,%d,%d]]\n",tet->vertices[2][0],tet->vertices[2][1],tet->vertices[2][2]);
}

/*
 * Returns whether the given triangle lies in a boundary plane of the cub
 * given by dimensions dim
 */
 
#define mat3_col_equal(mat,col,val) (mat[0][col] == val && mat[1][col] == val && mat[2][col] == val)
int triangle_boundary_cube(ptriangle triang, int dim) {
  return (mat3_col_equal(triang->vertices,0,0) || mat3_col_equal(triang->vertices,0,dim) ||
          mat3_col_equal(triang->vertices,1,0) || mat3_col_equal(triang->vertices,1,dim) ||
          mat3_col_equal(triang->vertices,2,0) || mat3_col_equal(triang->vertices,2,dim));
}

#define arr3_sum(arr) (arr[0] + arr[1] + arr[2])
int triangle_boundary_tet(ptriangle triang, int dim) {
  return (mat3_col_equal(triang->vertices,0,0) || //X = 0 plane
          mat3_col_equal(triang->vertices,1,0) || //Y = 0 plane
          mat3_col_equal(triang->vertices,2,0) || //Z = 0 plane
         (arr3_sum(triang->vertices[0]) == dim &&
          arr3_sum(triang->vertices[1]) == dim &&
          arr3_sum(triang->vertices[2]) == dim)); //X + Y + Z = dim plane
  
}
void triangle_symmetry(ptriangle triang, int sym,int dim, ptriangle res) {
  apply_symmetry(sym, dim, triang->vertices[0], res->vertices[0]);
  apply_symmetry(sym, dim, triang->vertices[1], res->vertices[1]);
  apply_symmetry(sym, dim, triang->vertices[2], res->vertices[2]);
}

triangle_list acute_triangle_recur(arr3 dim){
  /*triangle_list result = {NULL, 0, {dim[0],dim[1],dim[2]}};
  triangle_list prev;
  int maxdim, axis;
  size_t count = 0;
  ptriangle t_arr = NULL;
  size_t len = 0;
  arr3 *cube_pts;
  maxdim = maxArr3(dim, &axis);
  if (maxdim == 0 || (dim[0] + 1) * (dim[1] + 1) * (dim[2]+1) < 3)
    return result;
  arr3 new_dim = {dim[0],dim[1],dim[2]};
  new_dim[axis]--;
  prev = acute_triangle_recur(new_dim);
   
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
  //Loop through all combinations
  do {
    size_t i = arr[0], j = arr[1], k = arr[2];
    //Find all acute tetrahedra with two vertices on both "edge-planes"
    if ((cube_pts[i][axis] == 0 || cube_pts[j][axis] == 0 || 
       cube_pts[k][axis] == 0) &&
      (cube_pts[i][axis] == maxdim || cube_pts[j][axis] == maxdim ||
       cube_pts[k][axis] == maxdim))
    {
      triangle cur_tet = (triangle) {{{cube_pts[i][0],cube_pts[i][1],cube_pts[i][2]},
                                {cube_pts[j][0],cube_pts[j][1],cube_pts[j][2]},
                                {cube_pts[k][0],cube_pts[k][1],cube_pts[k][2]}}};
      if (mat3_triangle_acute(cur_tet.vertices)) {
        count++;
        t_arr = (ptriangle) realloc(t_arr, count * sizeof(triangle));
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
  ptriangle prev_edges = NULL;
  size_t edge_count = 0;
  for (int i = 0; i < prev.len; i++) {
    if (prev.t_arr[i].vertices[0][axis] == maxdim - 1 ||
        prev.t_arr[i].vertices[1][axis] == maxdim - 1 ||
        prev.t_arr[i].vertices[2][axis] == maxdim - 1) //Triangle was on expanding edge-plane
    {
      edge_count++;
      prev_edges = (ptriangle) realloc(prev_edges, edge_count * sizeof(triangle));
      if (prev_edges == NULL) {
        puts("Error allocating memory!!");
        exit(1);
      }
      prev_edges[edge_count - 1] = prev.t_arr[i];
      //Shift this edge tetrahedra
      prev_edges[edge_count - 1].vertices[0][axis] += 1;
      prev_edges[edge_count - 1].vertices[1][axis] += 1;
      prev_edges[edge_count - 1].vertices[2][axis] += 1;
    }
  }
  //Now the total set of acute tetrahedra constis of the old set + old_edge shifted + new edge calculated above  
  t_arr = realloc(t_arr,(count + edge_count + prev.len) * sizeof(triangle));
  
  memcpy(t_arr + count , prev_edges, edge_count * sizeof(triangle));
  memcpy(t_arr + (count + edge_count), prev.t_arr, prev.len * sizeof(triangle));
  
  result.t_arr = t_arr;
  result.len = (count + edge_count + prev.len);
  free(cube_pts);
  free(prev_edges);
  free(prev.t_arr);
  
  printf("Finished dimension: %d,%d,%d\n",dim[0],dim[1],dim[2]);
  printf("  Previous dimension: %d,%d,%d\n",new_dim[0],new_dim[1],new_dim[2]);
  printf("  Triangles added: %zu\n", (count + edge_count));
  printf("  Triangles prev:  %d\n", (prev.len));
  printf("  Amount of triangles: %zu\n", combination(len,3));
  printf("  Fraction sharp: %f\n", (float) result.len / combination(len,3));
  printf("  Total triangles: %d\n\n", result.len);
  return result;  
  */
}
