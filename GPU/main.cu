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

#define DIM 150
#define THREADS_BLOCK 512



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
  /*
   * Might as well copy the points once
   */
  cur_pt[0] = idx / ((dim+1) * (dim + 1));
  t = idx % ((dim+1) * (dim+1));
  cur_pt[1] = t / (dim + 1);
  cur_pt[2] = t % (dim + 1);
  
  subArr3(triang->vertices[1], triang->vertices[0], P[0]);
  subArr3(triang->vertices[2], triang->vertices[0], P[1]);
  subArr3(cur_pt, triang->vertices[0], P[2]);
  subArr3(triang->vertices[2], triang->vertices[1], P[3]);
  subArr3(cur_pt, triang->vertices[1], P[4]); 

  crossArr3(P[4],P[3], normals[0]); //Normal on facet 1,2,3  
  crossArr3(P[1],P[2], normals[1]); //Normal on facet 0,2,3
  crossArr3(P[2],P[0], normals[2]); //Normal on facet 0,1,3
  crossArr3(P[0],P[1], normals[3]); //Normal on facet 0,1,2
  //Normals[3] is the normal on the triangle plane
  t = dotArr3(normals[3], triang->vertices[0]); //Find the constant specific for this plane  
  
  result[idx] = ((dotArr3(normals[1], normals[2]) < 0) &
                 (dotArr3(normals[2], normals[3]) < 0) &
                 (dotArr3(normals[1], normals[3]) < 0) &
                 (dotArr3(normals[0], normals[1]) < 0) &
                 (dotArr3(normals[0], normals[2]) < 0) &
                 (dotArr3(normals[0], normals[3]) < 0)) //First bit
                 | 
                ((dotArr3(cur_pt, normals[3]) <  t) << 1); //second bit
                 

}


int facet_cube_acute_gpu(ptriangle triang, facet_acute_data * data, int mode, unsigned char * res_h) {
  /*
   * Every facet of an acute tetrahedron needs to be acute. If this facet is not even acute
   * we may directly stop checking this facet as it's never going to be part of an acute tetrahedron
   */
  if (!mat3_triangle_acute(triang->vertices)) 
    return 0;
    
    
  data->boundary_triangle = data->boundary_func(triang, data->cube->dim); //Boundary plane only needs acute tetra on 1 side
  int dim = data->cube->dim[0];
  data->acute_above = 0;
  data->acute_below = 0;
  data->tetra_above_len = 0;
  data->tetra_above = NULL;
  data->tetra_below_len = 0;
  data->tetra_below = NULL;
  
  
  
  for (size_t i = 0; i < data->cube->len; i++){
    if (res_h[i] == 1 && !data->acute_above ||
        res_h[i] == 3 && !data->acute_below) {
      arr3 cur_pt;
      int t;
      cur_pt[0] = i / ((dim+1) * (dim + 1));
      t = i % ((dim+1) * (dim+1));
      cur_pt[1] = t / (dim + 1);
      cur_pt[2] = t % (dim + 1);
      //All the facets must be in the acute_list
      if ((mode == FACET_ACUTE_LIST) && !facet_tetra_list(triang, cur_pt, data->acute_list))
        continue;   
      
      //Explicitly create a list of the acute tetrahedron
      if (mode == FACET_ACUTE_TETRA) {
        tetra test_tetra;
        memcpy(test_tetra.vertices, cur_pt, sizeof(arr3));
        memcpy(test_tetra.vertices + 1, triang->vertices, 3 * sizeof(arr3)); 
      
        if (res_h[i] == 1)
          tetra_add_array(test_tetra, &data->tetra_above, &data->tetra_above_len);
        else
          tetra_add_array(test_tetra, &data->tetra_below, &data->tetra_below_len);
      }
      //We only need to know if tetrahedron above and below acute
      else {
        if (res_h[i] == 1)
          data->acute_above = 1;
        else 
          data->acute_below = 1;
        if ((data->acute_above && data->acute_below) || data->boundary_triangle) {
          free(res_h);
          return 1;
        }
      }
    }
  }
  free(res_h);
  if (mode == FACET_ACUTE_TETRA) {
    data->acute_above = (data->tetra_above_len > 0);
    data->acute_below = (data->tetra_below_len > 0);
    if ((data->acute_above && data->acute_below) || (data->boundary_triangle && (data->acute_above || data->acute_below)))
      return 1;
  }
  return 0;  
}

int * facets_cube_acute_gpu(ptriangle triang, size_t n, facet_acute_data * data) {
  size_t len = data->cube->len;
  int * acute = (int *)  malloc(sizeof(int) * n);
  unsigned char * res_h, *res_d;
  ptriangle ptriang_d; 
  
  checkCudaCall(cudaMalloc(&res_d, len * sizeof(unsigned char)));
  checkCudaCall(cudaMalloc(&ptriang_d, sizeof(triangle))); 
  res_h = (unsigned char *) malloc(len * sizeof(unsigned char));  
  
  for (size_t i = 0; i < n; i++) {
    checkCudaCall(cudaMemcpy(ptriang_d, triang + i, sizeof(triangle), cudaMemcpyHostToDevice));
    tet_acute_kernel <<< len/THREADS_BLOCK + 1  , THREADS_BLOCK >>> (ptriang_d, data->cube->dim[0], res_d, len);
    checkCudaCall(cudaMemcpy(res_h, res_d, len * sizeof(unsigned char), cudaMemcpyDeviceToHost));
    acute[i] = facet_cube_acute_gpu(triang + i,data,FACET_ACUTE, res_h);
  }
  
  checkCudaCall(cudaFree(res_d));
  checkCudaCall(cudaFree(ptriang_d));
  return acute;

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
#define SIZE_LIST 100000
int main(void)
{
  clock_t begin, end;
  timeval t1,t2;
  triangle * triangles = (triangle * ) malloc(sizeof(triangle) * SIZE_LIST);
  int * acute;
  for (int i = 0; i < SIZE_LIST; i++)
    triangles[i] = rand_triangle(DIM);
  //triangle triang = rand_triangle(DIM);
  arr3 dim = {DIM,DIM,DIM};
  cube_points cube_pts = gen_cube_points(dim);
  facet_acute_data parameters;
  parameters.cube = &cube_pts;
  parameters.boundary_func = &triangle_boundary_cube;
  printf("Triangle: \n");
  gettimeofday(&t1,NULL);
  begin = clock();
  acute = facets_cube_acute_gpu(triangles,SIZE_LIST,&parameters);
  end = clock();
  printf("Acute_GPU: %d\n", acute);
  gettimeofday(&t2,NULL);
  printf("Wall time : %ld\n", ((t2.tv_sec * 1000000 + t2.tv_usec) - (t1.tv_sec * 1000000 + t1.tv_usec))/1000); 
  printf("Time taken on CPU: %f sec\n", float( (end - begin) )/ CLOCKS_PER_SEC);
  gettimeofday(&t1,NULL);
  begin = clock();
  for (int i =0; i < SIZE_LIST; i++)
    acute[i] = facet_cube_acute(triangles + i,&parameters,FACET_ACUTE);
  end = clock();
  printf("Acute: %d\n", acute);
  gettimeofday(&t2,NULL);
  printf("Wall time : %ld\n", ((t2.tv_sec * 1000000 + t2.tv_usec) - (t1.tv_sec * 1000000 + t1.tv_usec))/1000); 
  printf("Time taken on CPU: %f sec\n", float( (end - begin) )/ CLOCKS_PER_SEC);
  free(cube_pts.points);

}
