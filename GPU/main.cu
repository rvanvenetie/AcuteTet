// example1.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <iostream>
#include <cuda.h>
extern "C" {
#include "../vector.h"
#include "../triangle.h"  
#include "../tetraeder.h"
  
}
using namespace std;

/* Utility function, use to do error checking.

   Use this function like this:

   checkCudaCall(cudaMalloc((void **) &deviceRGB, imgS * sizeof(color_t)));

   And to check the result of a kernel invocation:

   checkCudaCall(cudaGetLastError());
*/
static void checkCudaCall(cudaError_t result) {
    if (result != cudaSuccess) {
        cerr << "cuda error: " << cudaGetErrorString(result) << endl;
        exit(1);
    }
}

// Kernel that executes on the CUDA device
__global__ void square_array(float *a, int N)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx<N) a[idx] = a[idx] * a[idx];
}

__global__ void tet_acute_kernel(triangle triang, arr3 * new_vec, int * result, size_t N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx>= N)
    return;
  
  arr3 P[5]; 
  arr3 normals[4]; 
  subArr3(triang.vertices[1], triang.vertices[0], P[0]);
  subArr3(triang.vertices[2], triang.vertices[0], P[1]);
  subArr3(triang.vertices[2], triang.vertices[1], P[3]);
  subArr3(new_vec[idx], triang.vertices[1], P[4]); 
  subArr3(new_vec[idx], triang.vertices[0], P[2]);
  
  crossArr3(P[2],P[0], normals[2]); //Normal on facet 0,1,3
  crossArr3(P[0],P[1], normals[3]); //Normal on facet 0,1,2
  crossArr3(P[1],P[2], normals[1]); //Normal on facet 0,2,3
  crossArr3(P[4],P[3], normals[0]); //Normal on facet 1,2,3

  result[idx] = (dotArr3(normals[1], normals[2]) < 0 &&
                 dotArr3(normals[1], normals[3]) < 0 &&
                 dotArr3(normals[0], normals[1]) < 0 &&
                 dotArr3(normals[0], normals[2]) < 0 &&
                 dotArr3(normals[0], normals[3]) < 0);
  
}



triangle rand_triangle(int dim) {
  triangle result;
  result.vertices[0][0] = rand() % dim;
  result.vertices[0][1] = rand() % dim;
  result.vertices[0][2] = rand() % dim;
  result.vertices[1][0] = rand() % dim;
  result.vertices[1][1] = rand() % dim;
  result.vertices[1][2] = rand() % dim;
  result.vertices[2][0] = rand() % dim;
  result.vertices[2][1] = rand() % dim;
  result.vertices[2][2] = rand() % dim;
  return result;
}

#define DIM 7

int main(void)
{
  triangle triang = rand_triangle(DIM);
  arr3 dim = {DIM,DIM,DIM};
  cube_points cube_pts = gen_cube_points(dim);
  
  arr3 * cube_d;
  int * res_h, *res_d;
  
  checkCudaCall(cudaMalloc(&cube_d, cube_pts.len * sizeof(arr3)));
  checkCudaCall(cudaMalloc(&res_d, cube_pts.len * sizeof(int)));
  
  res_h = (int *) malloc(cube_pts.len * sizeof(int));
  for (int i = 0; i < 5; i++)
    printf("%d %d\n", i, res_h[i]);
  
  checkCudaCall(cudaMemcpy(cube_d, &cube_pts.points[0], cube_pts.len * sizeof(arr3), cudaMemcpyHostToDevice));
  
  tet_acute_kernel <<< 1, cube_pts.len >>> (triang, cube_d, res_d, cube_pts.len);
  
  checkCudaCall(cudaGetLastError());
  checkCudaCall(cudaMemcpy(res_h, res_d, cube_pts.len * sizeof(int), cudaMemcpyDeviceToHost));

  for (int i=0; i<5; i++) {
    printf("Cuda: %d %d\n", i, res_h[i]);
    tetra test_tetra;
    memcpy(test_tetra.vertices, &cube_pts.points[i], sizeof(arr3));
    memcpy(test_tetra.vertices + 1, triang.vertices, 3 * sizeof(arr3)); 
    printf("Normal: %d %d\n", i , tetra_acute(&test_tetra));
  }
   
  cudaFree(cube_d);
  cudaFree(res_d);
  free(cube_pts.points);
  free(res_h);
}
