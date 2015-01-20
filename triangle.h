#ifndef TRIANGLE_H
#define TRIANGLE_H
#include <inttypes.h>
#include "vector.h"
//#include "datastructures.h"

typedef struct triangle
{
  arr3 vertices[3];
} triangle, *ptriangle;

typedef struct triangle_2d
{
  arr2 vertices[3];
} triangle_2d, *ptriangle_2d;

typedef struct triangle_list
{
  ptriangle t_arr;
  int len;
  int dim;   
} triangle_list;

//Returns area times two
int triangle_area_x2(arr2 v1, arr2 v2, arr2 v3);
int triangle_acute_2d(arr2 v1, arr2 v2, arr2 v3);

void arr3_to_triangle(arr3 v0, arr3 v1, arr3 v2, ptriangle triang);
int tri_tri_equal(ptriangle t1, ptriangle t2);
int arr3_triangle_acute(arr3 v0, arr3 v1, arr3 v2);
int triangle_acute(ptriangle triang);
int triangle_boundary_cube(ptriangle triang, int dim);
int triangle_boundary_tet(ptriangle triang, int dim);

void triangle_normal(ptriangle triang, arr3 normal);
void print_triangle(ptriangle tet);
void print_triangle_2d(ptriangle_2d tet);
void triangle_symmetry(ptriangle triang, int sym,int dim, ptriangle result);

/*
 * Calculates the edge vectors of the triangle given by 3 points: v0, v1, v2.
 */
#define triangle_sides(v0,v1,v2,tmp) { \
  subArr3(v1,v0, tmp[0]); \
  subArr3(v2,v0, tmp[1]); \
  subArr3(v2,v1, tmp[2]); }
  
/*
 * Given the three edge vectors of a triangle, this returns if the triangle is acute.
 * Convention is that P2 must point from the end of P0 to end of P1.
 */
#define triangle_sides_acute(P0,P1,P2) ((dotArr3(P0,P1) > 0) && \
          (dotArr3(P1,P2) > 0) && \
          (dotArr3(P0,P2) < 0))

/*
 * Does the same as above, but now P is an array of arr3
 */           
#define triangle_P_acute(P) (triangle_sides_acute(P[0],P[1],P[2]))


#endif
