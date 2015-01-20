#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vector.h"
#include "triangle.h"

//Triangle with vertices (edge->vertices[0], edge->vertices[1], apex)
int triangle_acute_2d(arr2 v1, arr2 v2, arr2 v3) {
	arr2 P[3]; //Edges
	subArr2(P[0], v2, v1); //b - a
	subArr2(P[1], v3, v1); //c - a
	subArr2(P[2], v3, v2); //c  -b

	return ((dotArr2(P[0],P[1]) > 0) &&
			    (dotArr2(P[1],P[2]) > 0) &&
		      (dotArr2(P[0],P[2]) < 0));
}
int triangle_area_x2(arr2 vA, arr2 vB, arr2 vC) {
	int //Easy of notation
		XA = vA[0], XB = vB[0], XC = vC[0],
		YA = vA[1], YB = vB[1], YC = vC[1];

	return abs((XB - XA)*(YC - YA) - (XC-XA) * (YB - YA));
}

void arr3_to_triangle(arr3 v0, arr3 v1, arr3 v2, ptriangle triang) {
  copyArr3(triang->vertices[0], v0);
  copyArr3(triang->vertices[1], v1);
  copyArr3(triang->vertices[2], v2);
}

int tri_tri_equal(ptriangle t1, ptriangle t2) {
  return (vert_vert_share_count(t1->vertices, 3, t2->vertices, 3) == 3);
}

int arr3_triangle_acute(arr3 v0, arr3 v1, arr3 v2) {
  arr3 P[3];
  triangle_sides(v0,v1,v2,P);
  return triangle_P_acute(P);
}

int triangle_acute(ptriangle triang) {
   return arr3_triangle_acute(triang->vertices[0], triang->vertices[1], triang->vertices[2]);  
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

void print_triangle_2d(ptriangle_2d tet) {
  printf("[[%d,%d],\n",tet->vertices[0][0],tet->vertices[0][1]);
  printf(" [%d,%d],\n",tet->vertices[1][0],tet->vertices[1][1]);
  printf(" [%d,%d]]\n",tet->vertices[2][0],tet->vertices[2][1]);
}
/*
 * Returns whether the given triangle lies in a boundary plane of the cube
 * given by dimensions dim.
 */
 
#define mat3_col_equal(mat,col,val) (mat[0][col] == val && mat[1][col] == val && mat[2][col] == val)
int triangle_boundary_cube(ptriangle triang, int dim) {
  return (mat3_col_equal(triang->vertices,0,0) || mat3_col_equal(triang->vertices,0,dim) ||
          mat3_col_equal(triang->vertices,1,0) || mat3_col_equal(triang->vertices,1,dim) ||
          mat3_col_equal(triang->vertices,2,0) || mat3_col_equal(triang->vertices,2,dim));
}

/*
 * Returns whether the given triangle lies on the boundary of the unit tetrahedron
 * of dimensions dim.
 */
#define arr3_sum(arr) (arr[0] + arr[1] + arr[2])
int triangle_boundary_tet(ptriangle triang, int dim) {
  return (mat3_col_equal(triang->vertices,0,0) || //X = 0 plane
          mat3_col_equal(triang->vertices,1,0) || //Y = 0 plane
          mat3_col_equal(triang->vertices,2,0) || //Z = 0 plane
         (arr3_sum(triang->vertices[0]) == dim &&
          arr3_sum(triang->vertices[1]) == dim &&
          arr3_sum(triang->vertices[2]) == dim)); //X + Y + Z = dim plane
  
}

/*
 * Applies symmetry on each vertex of the triangle
 */
void triangle_symmetry(ptriangle triang, int sym,int dim, ptriangle res) {
  apply_symmetry(sym, dim, triang->vertices[0], res->vertices[0]);
  apply_symmetry(sym, dim, triang->vertices[1], res->vertices[1]);
  apply_symmetry(sym, dim, triang->vertices[2], res->vertices[2]);
}

