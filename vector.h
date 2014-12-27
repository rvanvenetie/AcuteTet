#ifndef VECTOR_H
#define VECTOR_H

#include <stdlib.h>
typedef struct vec3
{
     int pt[3];
} vec3,*pvec3;

typedef  int arr2[2];
typedef  int arr3[3];

typedef struct cube_points
{
  arr3 * points;
  int dim;
  size_t len;
} cube_points;


typedef unsigned short vert_index;
typedef unsigned short  ***vert_index_array;

#ifdef INLINE_MACROS
  #define subArr3(u,v, result) {\
    result[0] = u[0] - v[0]; \
    result[1] = u[1] - v[1]; \
    result[2] = u[2] - v[2];}

  #define crossArr3(u,v, result) {\
    result[0] = u[1] * v[2] - u[2] * v[1]; \
    result[1] = u[2] * v[0] - u[0] * v[2]; \
    result[2] = u[0] * v[1] - u[1] * v[0];}
 


  #define dotArr3(u,v) (u[0]*v[0] + u[1]*v[1] + u[2]*v[2])
  
  #define zeroArr3(u) (u[0] ==0 && u[1] == 0 && u[2] == 0)
  #define equalArr3(u,v) (u[0] == v[0] && u[1] == v[1] && u[2] == v[2])
 
  #define mirror_sym {t = x;x = y;y = t;}    
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
  #define apply_symmetry(sym,dim,pt, result) {\
    int t; \
    int x = pt[0], y = pt[1], z = pt[2]; \
    switch (sym) { \
      case 0:  face0; rot0;   break;     \
      case 1:  face0; rot90;  break;     \
      case 2:  face0; rot180; break;     \
      case 3:  face0; rot270; break;     \
      case 4:  face1; rot0;   break;     \
      case 5:  face1; rot90;  break;     \
      case 6:  face1; rot180; break;     \
      case 7:  face1; rot270; break;     \
      case 8:  face2; rot0;   break;     \
      case 9:  face2; rot90;  break;     \
      case 10: face2; rot180; break;     \
      case 11: face2; rot270; break;     \
      case 12: face3; rot0;   break;     \
      case 13: face3; rot90;  break;     \
      case 14: face3; rot180; break;     \
      case 15: face3; rot270; break;     \
      case 16: face4; rot0;   break;     \
      case 17: face4; rot90;  break;     \
      case 18: face4; rot180; break;     \
      case 19: face4; rot270; break;     \
      case 20: face5; rot0;   break;     \
      case 21: face5; rot90;  break;     \
      case 22: face5; rot180; break;     \
      case 23: face5; rot270; break;     \
      case 24: mirror_sym; face0; rot0;   break;\
      case 25: mirror_sym; face0; rot90;  break;\
      case 26: mirror_sym; face0; rot180; break;\
      case 27: mirror_sym; face0; rot270; break;\
      case 28: mirror_sym; face1; rot0;   break;\
      case 29: mirror_sym; face1; rot90;  break;\
      case 30: mirror_sym; face1; rot180; break;\
      case 31: mirror_sym; face1; rot270; break;\
      case 32: mirror_sym; face2; rot0;   break;\
      case 33: mirror_sym; face2; rot90;  break;\
      case 34: mirror_sym; face2; rot180; break;\
      case 35: mirror_sym; face2; rot270; break;\
      case 36: mirror_sym; face3; rot0;   break;\
      case 37: mirror_sym; face3; rot90;  break;\
      case 38: mirror_sym; face3; rot180; break;\
      case 39: mirror_sym; face3; rot270; break;\
      case 40: mirror_sym; face4; rot0;   break;\
      case 41: mirror_sym; face4; rot90;  break;\
      case 42: mirror_sym; face4; rot180; break;\
      case 43: mirror_sym; face4; rot270; break;\
      case 44: mirror_sym; face5; rot0;   break;\
      case 45: mirror_sym; face5; rot90;  break;\
      case 46: mirror_sym; face5; rot180; break;\
      case 47: mirror_sym; face5; rot270; break;\
    }\
    result[0] = x;result[1] = y;result[2] = z;}  
    
#else
  void crossArr3(arr3 u, arr3 v, arr3 result);
  void subArr3(arr3 u, arr3 v, arr3 result);
  int dotArr3(arr3 u, arr3 v);
  void apply_symmetry(int sym,int dim,arr3 pt, arr3 result);
  int zeroArr3(arr3 u);
  int equalArr3(arr3 u,arr3 v); 
#endif

int vert_vert_share_count(arr3 * vert1, int len1, arr3 * vert2, int len2);
int vert_in_plane_count(arr3 normal, int d, arr3 * pts, int len_pts);
void negArr3(arr3 result);
void printVector(vec3 u);
int dotVector(vec3 u, vec3 v);
vec3 negVector(vec3 u);
vec3 scalarVector(int alpha, vec3 u);
vec3 crossVector(vec3 u, vec3 v);
vec3 addVectors(vec3 u, vec3 v);
vec3 subVectors(vec3 u, vec3 v);
int maxArr3(arr3 u, int * axis);
void copyArr3(arr3 dest, arr3 source);
void printArr3(arr3 u);

cube_points gen_tet_points(int dim);
cube_points gen_cube_points(int dim);
cube_points gen_fund_points(int dim);
cube_points gen_cube_sparse_points(int dim);
cube_points gen_fund_sparse_points(int dim);
void gen_sparse_axis(int dim, int ** axis, size_t * len);

void randomArr3(int dim, arr3 result);
#endif
