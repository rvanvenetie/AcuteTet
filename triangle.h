#ifndef TRIANGLE_H
#define TRIANGLE_H
#include <inttypes.h>
#include "vector.h"
//#include "datastructures.h"

typedef struct Triangle
{
    vec3 vertices[3];
    vec3 P[3];
} Triangle;

typedef struct triangle
{
  arr3 vertices[3];
} triangle, *ptriangle;

typedef struct triangle_list
{
  ptriangle t_arr;
  int len;
  arr3 dim;   
} triangle_list;


Triangle triange_init(vec3 v1, vec3 v2, vec3 v3);


int triangle_acute(Triangle *triang);
int mat3_triangle_acute(mat3 v);
int arr3_triangle_acute(arr3 v0, arr3 v1, arr3 v2);
int triangle_boundary_cube(ptriangle triang, int dim);
int triangle_boundary_tet(ptriangle triang, int dim);

void triangle_normal(ptriangle triang, arr3 normal);
void print_triangle(ptriangle tet);
triangle_list acute_triangle_recur(arr3 dim);
void triangle_symmetry(ptriangle triang, int sym,int dim, ptriangle result);


//Puts the sides of this triangle in the variable tmp.
#define triangle_sides(v0,v1,v2,tmp) { \
  subArr3(v1,v0, tmp[0]); \
  subArr3(v2,v0, tmp[1]); \
  subArr3(v2,v1, tmp[2]); }
  
/*
 * Given the three sides of a triangle, this returns if the triangle is acute.
 * Convention is that P2 must point from the end of P0 to end of P1.
 */
#define triangle_sides_acute(P0,P1,P2) ((dotArr3(P0,P1) > 0) && \
          (dotArr3(P1,P2) > 0) && \
          (dotArr3(P0,P2) < 0))
           
#define triangle_P_acute(P) (triangle_sides_acute(P[0],P[1],P[2]))
#endif
