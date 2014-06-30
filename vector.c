/* Computer Graphics, Assignment, Intersections
 * Filename ........ normals.c
 * Description ..... Calculates Plane/Vertex normals
 * Date ............ 19.08.2008
 *
 * Student name Ruben Brokkelkamp & Raymond van Venetië
 * Student email rubenbrokkelkamp@gmail.com & rayvanve@hotmail.com
 * Collegekaart 6282628 & 6303900
 * Date 15-09-2011
 * Comments ........
 *
 * (always fill in these fields before submitting!!)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vector.h"


#ifndef INLINE_MACROS
void subArr3(arr3 u, arr3 v, arr3 result) {
  result[0] = u[0] - v[0];
  result[1] = u[1] - v[1];
  result[2] = u[2] - v[2];
}

int dotArr3(arr3 u, arr3 v){
  return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];  
}
void crossArr3(arr3 u, arr3 v, arr3 result){
  result[0] = u[1] * v[2] - u[2] * v[1];
  result[1] = u[2] * v[0] - u[0] * v[2];
  result[2] = u[0] * v[1] - u[1] * v[0];
}
int zeroArr3(arr3 u) {
  return (u[0] == 0 && u[1] == 0 && u[2] == 0);
}

int equalArr3(arr3 u,arr3 v){
  return (u[0] == v[0] && u[1] == v[1] && u[2] == v[2]);
}
#endif

void negArr3(arr3 result) {
  result[0] = -result[0];
  result[1] = -result[1];
  result[2] = -result[2];
}

void printVector(vec3 u){
  printf("[%d,%d,%d]\n",u.pt[0],u.pt[1],u.pt[2]);
}
vec3 scalarVector(int alpha, vec3 u){
  return ((vec3) {{alpha * u.pt[0], alpha * u.pt[1], alpha * u.pt[2]}});
}

vec3 negVector(vec3 u){
  return ((vec3) {{-u.pt[0],-u.pt[1],-u.pt[2]}});
}

int dotVector(vec3 u, vec3 v){
  return u.pt[0] * v.pt[0] + u.pt[1] * v.pt[1] + u.pt[2]*v.pt[2];
}


vec3 crossVector(vec3 u, vec3 v){
  return ((vec3) 
    {{u.pt[1] * v.pt[2] - u.pt[2] * v.pt[1],
      u.pt[2] * v.pt[0] - u.pt[0] * v.pt[2],
      u.pt[0] * v.pt[1] - u.pt[1] * v.pt[0]}});
}
vec3 addVectors(vec3 u, vec3 v){
  return ((vec3) 
  {{ u.pt[0] + v.pt[0],
     u.pt[1] + v.pt[1],
     u.pt[2] + v.pt[2]}});
}

vec3 subVectors(vec3 u, vec3 v){
  vec3 result;
  result.pt[0] = u.pt[0] - v.pt[0];
  result.pt[1] = u.pt[1] - v.pt[1];
  result.pt[2] = u.pt[2] - v.pt[2];
  return result;
  /*
  return ((vec3) 
  {{ u.pt[0] - v.pt[0],
     u.pt[1] - v.pt[1],
     u.pt[2] - v.pt[2]}});
  */
}



int maxArr3(arr3 u, int * axis) {
  if (u[0] >= u[1] && u[0] >= u[2]) {
    *axis = 0;
    return u[0];
  } else if(u[1] >= u[0] && u[1] >= u[2]) {
    *axis = 1;
    return u[1];
  } else {
    *axis = 2;
    return u[2];
  }
}

void printArr3(arr3 u) {
  printf("[%d,%d,%d]\n",u[0],u[1],u[2]);   
}
void initMat3(vec3 c1, vec3 c2, vec3 c3,mat3 result){
  memcpy(result[0],c1.pt,sizeof(c1.pt));
  memcpy(result[1],c2.pt,sizeof(c1.pt));
  memcpy(result[2],c3.pt,sizeof(c1.pt));
}

void transMat3(mat3 mat, mat3 result) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      result[i][j] = mat[j][i];
}

vec3 mat3Col(mat3 mat, int col) {
  return ((vec3) {{mat[0][col], 
                   mat[1][col],
                   mat[2][col]}});
}

vec3 mat3Row(mat3 mat, int row) {
  return ((vec3) {{mat[row][0], mat[row][1], mat[row][2]}});
}
/*
 * Matrix times vector
 */
vec3 mat3Vector(mat3 mat, vec3 v) {
  vec3 result = (vec3) {{0}};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      result.pt[i] += mat[i][j] * v.pt[j];
  return result;
}
/*
 * Calculates MN as matrices
 */
void dotMat3(mat3 M,mat3 N, mat3 result) {
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++){
      result[i][j]=0;
      for(int k = 0; k < 3; k++)
          result[i][j] += M[i][k]*N[k][j];
    }
}

void printMatrix(mat3 mat){
  printf("[[%d,%d,%d],\n",mat[0][0],mat[0][1],mat[0][2]);
  printf(" [%d,%d,%d],\n",mat[1][0],mat[1][1],mat[1][2]);
  printf(" [%d,%d,%d]]\n",mat[2][0],mat[2][1],mat[2][2]);
}


/*
 * De symmetrie group is gelijk aan S_4 x Z_2 (waarbij Z_2 de reflecties voorstellen en S_4 de rotaties)
 * We kunnen hem verder onderverdelen met S_4 = V_4 x 
 * 
 */
 
void symmetries_old(int sym,int dim,arr3 pt, arr3 result) {
  int t,face,rot,mirror,symmetry = sym;
  int x = pt[0], y = pt[1], z = pt[2];
  //Numbers above 24 mean mirroring
  mirror = symmetry / 24;
  //New symmetry is modulus 24
  symmetry = symmetry % 24;
  //Now symmetry 0..5, 6..11, 12..17, 18..24 give the face and finally rotation
  face = symmetry % 6;
  rot =  symmetry / 6;
  //printf("Symmetry %d, mir: %d, face: %d, rot: %d\n", sym, mirror, face, rot);
  /*
   * Mirror in XY plane
   */
  if (mirror) {
    t = x;
    x = y;
    y = t;    
  }
  /*
   * Rotate face 'down' to the bottom
   */
  switch (face) {
    case 0: //Z = 0 face, do nothing 
      break;
    case 1: //X = 0 face
      /*
       * X-coordinate stays the same
       * Y -> dim - z
       * Z -> Y
       */
      t = y;
      y = dim - z;
      z = t;
      break;
    case 2: //Y = 0 face
      t = x;
      x = dim - z;
      z = t;
      break;
    case 3: //Z = dim face, rotate over y-axis
      //X stays the same
      y = dim - y;
      z = dim - z;
      break;
    case 4: //X = dim face
      //X stays the same
      t = y;
      y = z;
      z = dim - t;
      break;
    case 5: //Y = dim face
      t = x;
      x = z;
      z = dim - t;
      break;
  }
  /*
   * Rotate with axis from top-to-bottom, 4 options, does not change z-axis
   */
  switch (rot) {
    case 0: //Do nothing
      break;
    case 1: //Rotate 90 degrees
      t = x;
      x = y;
      y = dim - t;
      break;
    case 2: //Rotate 180 degrees
      x = dim - x;
      y = dim - y;
      break;
    case 3: //Rotate 270 degrees
      t = x;
      x = dim - y;
      y = t;
      break;
  }
  result[0] = x;
  result[1] = y;
  result[2] = z;
}

  
#define mirror {t = x;x = y;y = t;}    
#define face0 {}
#define face1 {t = y;y = dim - z;z = t;}
#define face2 {t = x;x = dim - z;z = t;}
#define face3 {y = dim - y;z = dim - z;}
#define face4 {t = y;y = z;z = dim - t;}
#define face5 {t = x;x = z;z = dim - t;}
#define rot0 {}
#define rot90 {t = x; x = y; y = dim - t;}
#define rot180 {x = dim - x;y = dim - y;}
#define rot270 {t = x;x = dim - y;y = t;}

#ifndef INLINE_MACROS
  void apply_symmetry(int sym,int dim,arr3 pt, arr3 result) {
    int t;
    int x = pt[0], y = pt[1], z = pt[2];
    switch (sym) {
      case 0:  face0; rot0;   break;
      case 1:  face0; rot90;  break;
      case 2:  face0; rot180; break;
      case 3:  face0; rot270; break;
      case 4:  face1; rot0;   break;
      case 5:  face1; rot90;  break;
      case 6:  face1; rot180; break;
      case 7:  face1; rot270; break;
      case 8:  face2; rot0;   break;
      case 9:  face2; rot90;  break;
      case 10: face2; rot180; break;
      case 11: face2; rot270; break;
      case 12: face3; rot0;   break;
      case 13: face3; rot90;  break;
      case 14: face3; rot180; break;
      case 15: face3; rot270; break;
      case 16: face4; rot0;   break;
      case 17: face4; rot90;  break;
      case 18: face4; rot180; break;
      case 19: face4; rot270; break;
      case 20: face5; rot0;   break;
      case 21: face5; rot90;  break;
      case 22: face5; rot180; break;
      case 23: face5; rot270; break;
      case 24: mirror; face0; rot0;   break;
      case 25: mirror; face0; rot90;  break;
      case 26: mirror; face0; rot180; break;
      case 27: mirror; face0; rot270; break;
      case 28: mirror; face1; rot0;   break;
      case 29: mirror; face1; rot90;  break;
      case 30: mirror; face1; rot180; break;
      case 31: mirror; face1; rot270; break;
      case 32: mirror; face2; rot0;   break;
      case 33: mirror; face2; rot90;  break;
      case 34: mirror; face2; rot180; break;
      case 35: mirror; face2; rot270; break;
      case 36: mirror; face3; rot0;   break;
      case 37: mirror; face3; rot90;  break;
      case 38: mirror; face3; rot180; break;
      case 39: mirror; face3; rot270; break;
      case 40: mirror; face4; rot0;   break;
      case 41: mirror; face4; rot90;  break;
      case 42: mirror; face4; rot180; break;
      case 43: mirror; face4; rot270; break;
      case 44: mirror; face5; rot0;   break;
      case 45: mirror; face5; rot90;  break;
      case 46: mirror; face5; rot180; break;
      case 47: mirror; face5; rot270; break;         
    }
    result[0] = x;result[1] = y;result[2] = z;
  }
#endif  


/*
 * Generates all points inside the unit tetrahedra.
 * Total amount of points is given by:
 *  - For triangle of dim k = \sum_{n=0}^k (k - n + 1) = 1/2(1+k)(2+k)
 *  - For all triangles = \sum_{k=0}^dim 1/2(1+k)(2+k) = 1/6 (n+1)(n+2)(n+3)
 * vertex_index will hold an array used to reverse a point to it's index.
 * vertex_index[x][y][z] will give it's index in the points array.
 */
cube_points gen_tet_points(int dim) {
  cube_points result = {NULL, 
                        {dim,dim,dim},
                        (dim+1)*(dim+2)*(dim+3)/6};
  result.points = malloc(result.len *sizeof(arr3));
  int c = 0;
  
  for (int x = 0; x <= dim; x++)
    for (int y = 0; y <= dim - x; y++) 
      for (int z = 0; z <= dim - x -y; z++) {
        result.points[c][0] = x;
        result.points[c][1] = y;
        result.points[c][2] = z;
        c++;
      }

  return result;
}

cube_points gen_cube_points(int dim) {
  cube_points result = {NULL, 
                        dim,
                        ((dim + 1) * (dim + 1) * (dim+1))};
  result.points = malloc(result.len *sizeof(int) * 3);
  int c = 0;
  for (int x = 0; x <= dim; x++)
    for (int y = 0; y <= dim; y++)
      for (int z = 0; z <= dim; z++) {
        result.points[c][0] = x;
        result.points[c][1] = y;
        result.points[c][2] = z;
        c++;
      }
  return result;
}


/*
 * Let M = dim/2, the middle of an edge.
 * Take O = (M, M, M) the middle of the cube.
 * Take A = (0,0,0), B = (M, 0, 0) and C = (M,M,0)
 * We take the subset of points inside the tetrahedron ABC O.
 * Area of base = 1/2 * M * M
 * Volume of tetrahedron = 1/3 * M * (1/2 * M * M) = 1/3 * 1/2 * (1/2)^3 * dim^3 
 *  = 1/48 * dim^3 
 */
cube_points gen_fund_points(int dim){
  cube_points result = {NULL, 
                        {dim, dim, dim},
                        0};

  /*
   * Possible just loop through all points and check if inside?
   */
  int M = dim / 2; //dim has points 0..dim, automatically floored. As point 3.5 is inside
  for (int z = 0; z <= M; z++)
    //We now have triangle (z,z,z) - (z,m,z) - (m,m,z)
    for (int x = z; x <= M; x++) //Loop over x-axis
      for (int y = z; y <= x; y++) {
        result.points = realloc(result.points, (result.len+1) * sizeof(arr3));
        result.points[result.len][0] = x;
        result.points[result.len][1] = y;
        result.points[result.len][2] = z;
        result.len++;
      }
  
  return result;
}

void randomArr3(int dim, arr3 result) {
  result[0] = rand() % dim;
  result[1] = rand() % dim;
  result[2] = rand() % dim;
}

void copyArr3(arr3 dest, arr3 source) {
  dest[0] = source[0];
  dest[1] = source[1];
  dest[2] = source[2];
  
}
