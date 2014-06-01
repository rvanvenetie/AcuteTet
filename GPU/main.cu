// example1.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <time.h>
#include <sys/time.h>

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

__global__ void tet_acute_kernel(ptriangle triang, int dim, unsigned char * result, size_t N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx>= N)
    return;
  arr3 P[5]; 
  arr3 cur_pt;
  arr3 normals[4];
  int t; //Temp
  cur_pt[0] = idx / ((dim+1) * (dim + 1));
  t = idx % ((dim+1) * (dim+1));
  cur_pt[1] = t / (dim + 1);
  cur_pt[2] = t % (dim + 1);
  
  subArr3(triang->vertices[1], triang->vertices[0], P[0]);
  subArr3(triang->vertices[2], triang->vertices[0], P[1]);
  subArr3(cur_pt, triang->vertices[0], P[2]);
  subArr3(triang->vertices[2], triang->vertices[1], P[3]);
  subArr3(cur_pt, triang->vertices[1], P[4]); 
  
  crossArr3(P[2],P[0], normals[2]); //Normal on facet 0,1,3
  crossArr3(P[0],P[1], normals[3]); //Normal on facet 0,1,2
  crossArr3(P[1],P[2], normals[1]); //Normal on facet 0,2,3
  crossArr3(P[4],P[3], normals[0]); //Normal on facet 1,2,3
  result[idx] = (dotArr3(normals[1], normals[2]) < 0 &
                 dotArr3(normals[2], normals[3]) < 0 &
                 dotArr3(normals[1], normals[3]) < 0 &
                 dotArr3(normals[0], normals[1]) < 0 &
                 dotArr3(normals[0], normals[2]) < 0 &
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

#define DIM 150
#define THREADS_BLOCK 512
int main(void)
{
  clock_t begin, end;
  cudaEvent_t start, stop;
  float time_ms;
  triangle triang = rand_triangle(DIM);
  triangle * ptriang;
  arr3 dim = {DIM,DIM,DIM};
  cube_points cube_pts = gen_cube_points(dim);
  arr3 * cube_d;
  unsigned char * res_h, *res_d;
  timeval t1,t2;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  gettimeofday(&t1,NULL);
  begin = clock();
  checkCudaCall(cudaMalloc(&res_d, cube_pts.len * sizeof(unsigned char)));
  checkCudaCall(cudaMalloc(&ptriang, sizeof(triangle))); 
  res_h = (unsigned char *) malloc(cube_pts.len * sizeof(unsigned char));
  checkCudaCall(cudaMemcpy(ptriang, &triang, sizeof(triangle), cudaMemcpyHostToDevice));
  cudaEventRecord(start,0);
  tet_acute_kernel <<< cube_pts.len/THREADS_BLOCK + 1  , THREADS_BLOCK >>> (ptriang, DIM, res_d, cube_pts.len);
  cudaEventRecord(stop,0);  
  checkCudaCall(cudaGetLastError());
  checkCudaCall(cudaMemcpy(res_h, res_d, cube_pts.len * sizeof(unsigned char), cudaMemcpyDeviceToHost));
  checkCudaCall(cudaFree(res_d));
  end = clock();
  cudaEventElapsedTime(&time_ms,start,stop);
  gettimeofday(&t2,NULL);
  printf("Wall time : %ld\n", ((t2.tv_sec * 1000000 + t2.tv_usec) - (t1.tv_sec * 1000000 + t1.tv_usec))); 
  printf("Time taken on CPU: %f sec\n", float( (end - begin) )/ CLOCKS_PER_SEC);
  printf("Time taken op the GPU: %f msec\n", time_ms);
  gettimeofday(&t1,NULL);  
  begin = clock();
  for (int i=0; i<cube_pts.len; i++) {
    //printf("Cuda: %d %d\n", i, res_h[i]);
    tetra test_tetra;
    memcpy(test_tetra.vertices + 3, cube_pts.points + i, sizeof(arr3));
    memcpy(test_tetra.vertices , triang.vertices, 3 * sizeof(arr3)); 
    unsigned char  acute = (unsigned char) tetra_acute(&test_tetra);
    if (res_h[i] != acute) {
      printf("FAIL %d %d %d\n",i, res_h[i], acute );
    }
  }
  end  = clock();
  gettimeofday(&t2,NULL);
  printf("Time taken on the CPU: %f sec\n", float((end - begin)) / CLOCKS_PER_SEC);
  printf("Wall time : %ld\n", ((t2.tv_sec * 1000000 + t2.tv_usec) - (t1.tv_sec * 1000000 + t1.tv_usec))); 
  free(cube_pts.points);
  free(res_h);
}
