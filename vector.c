#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vector.h"


/*
 * Returns how many points of vert1 are also in vert2. If one of the arrays
 * contains multiple instances of the same point, all those instances are counted.
 */
int vert_vert_share_count(arr3 * vert1, int len1, arr3 * vert2, int len2) {
  int cnt = 0;
  for (int i  = 0; i < len1; i++)
    for (int j = 0; j < len2; j++)
      if (equalArr3(vert1[i], vert2[j]))
	cnt++;
  return cnt;
}
/*
 * Returns how many points are contained in the plane given by the normal and a point
 * in the plane.
 *
 * (Points in the plane satisfy equation dot(normal, pt) = dot(normal, pt_in_plane))
 * Where d = dot(normal, pt_in_plane)
 */
int vert_in_plane_count(arr3 normal, int d, arr3 * pts, int len_pts) {
  int cnt = 0;
  for (int i = 0; i < len_pts; i++)
    if (dotArr3(normal, pts[i]) == d)
      cnt++;

  return cnt;
}
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

void printArr2(arr2 u) {
  printf("[%d,%d]\n",u[0],u[1]);
}
void printArr3(arr3 u) {
  printf("[%d,%d,%d]\n",u[0],u[1],u[2]);   
}


/*
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
                        dim,
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
                        dim,
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

/*
 * Generates a sparse axis.
 */
void gen_sparse_axis(int dim, int ** axis, size_t * len) {
  (*axis) = calloc(dim + 1, sizeof(int));
  /*
   * We will first create an indicator array indicating what points are on the (*axis)
   * We will add the middle of the subintervals to the (*axis) (rounded towards the middle).
   * We only save the middle of the subinterval closest to the middle.
   */
  int left, right;
  (*axis)[dim / 2] = 1;
  //Start left of dim/2
  left = 0; right = dim / 2;
  while (left != right) {
    (*axis)[left] = 1;
    left = (left + right + 1) / 2; //Round up
  }
  //Start right of dim/2
  left = dim / 2; right = dim;
  while (left != right) {
    (*axis)[right] = 1;
    right = (left + right) / 2; //Round down
  }

  /*
   * Len will hold the amount of points in total
   */
  *len = 0;
  for (int i = 0; i <= dim; i++)
    if ((*axis)[i])
      (*len)++;
}
/*
 * Generates the points inside the cube on a `sparse' grid. This grid
 * is chosen in such a way that the amount of points in the centre is dense, compared to the
 * outher points.
 */
cube_points gen_cube_sparse_points(int dim) {
  int * axis;
  size_t len_axis;

  gen_sparse_axis(dim, &axis, &len_axis);
  cube_points result = {NULL, 
                        dim,
			(len_axis) * (len_axis) * (len_axis)};
  result.points = malloc(result.len *sizeof(arr3));

  int c = 0;
  for (int x = 0; x <= dim; x++) {
    if (!axis[x])
      continue;
    for (int y = 0; y <= dim; y++) {
      if (!axis[y])
        continue;
      for (int z = 0; z <= dim; z++)
      {
        if (!axis[z])
          continue;
        result.points[c][0] = x;
        result.points[c][1] = y;
        result.points[c][2] = z;
        c++;
      }
    }
  }
  free(axis);
  return result;
}

/*
 * Generates the fundamental points inside the cube with a sparse grid, see 
 * gen_fund_points
 */
cube_points gen_fund_sparse_points(int dim){
  cube_points result = {NULL, 
                        dim,
                        0};
  int * axis = NULL;
  size_t len_axis;
  gen_sparse_axis(dim, &axis, &len_axis);

  int M = dim / 2; //dim has points 0..dim, automatically floored. As point 3.5 is inside
  for (int z = 0; z <= M; z++) {
    if (!axis[z])
      continue;
    //We now have triangle (z,z,z) - (z,m,z) - (m,m,z)
    for (int x = z; x <= M; x++){ //Loop over x-axis
      if (!axis[x])
        continue;
      for (int y = z; y <= x; y++) {
        if (!axis[y])
          continue;

        result.points = realloc(result.points, (result.len+1) * sizeof(arr3));
        result.points[result.len][0] = x;
        result.points[result.len][1] = y;
        result.points[result.len][2] = z;
        result.len++;
      }
    }
  }
  free(axis);
  return result;
}
square_points gen_fund_square(int p) {
  square_points result = {NULL, 
                        p,
                        0};
	int m = p / 2;
	for (int x = 0; x <= m; x++)
		for (int y = 0; y <= x; y++)
		{
			result.points = realloc(result.points, (++result.len) * sizeof(arr2));
			result.points[result.len - 1][0] = x;
			result.points[result.len - 1][1] = y;
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
void copyArr2(arr2 dest, arr2 source) {
  dest[0] = source[0];
  dest[1] = source[1];
}
