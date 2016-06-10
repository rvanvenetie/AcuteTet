#pragma once
#include "vector.h"

// A memory set containing all the triangles
class TriangleSet {
  public:
    byte  *** t_arr; //Actual data array

    TriangleSet(int dim)
    {
    }

    // direct access to the vertices
    inline bool operator()(int a, int b, int c) const { 
      return t_arr[a][b][c/ 8] & (1 << (c % 8));
    }

};
// A memory set containing triangles, stored using the fundamental domain
class FundTriangleSet {
  public:
    FundTriangleSet(int dim) 
    {
    }

};

/*
tri_mem_list mem_list_fund_init(int dim, int init_value, int mode) {
  tri_mem_list result;
  memset(&result,0,sizeof(tri_mem_list));
  result.mode = mode;
  result.dim = dim; 

  cube_points cube, fund;
  if (mode == MEM_LIST_FUND)
  {
    gen_axis(dim, 0);
    cube = gen_cube_points(dim);
    fund = gen_fund_points(dim);
  } else {
    gen_axis(dim, 1); //HACKY
    cube = gen_cube_sparse_points(dim);
    fund = gen_fund_sparse_points(dim); 
  }
  result.mem_fund.cube_len = cube.len;
  result.mem_fund.fund_len = fund.len;
  //Generate vertex to index array
  gen_vertex_to_index_fund(axis_sparse_to_normal, axis_sparse_len, cube, fund, &result.mem_fund.vert_to_index);
  //Generate index array to vertex
  gen_vertex_from_index_fund(&result.mem_fund.vert_from_index, result.mem_fund.vert_to_index, cube);


  /* Initalize the 3D bit array.
   * First dimension only holds points in fundamental domain (the first fund.len points)
   * Second dimension then holds all the points bigger than first dimension (also outside fund domain)
   * Third dimension holds all points bigger than the third dimension
  unsigned char *** t_arr = calloc(fund.len, sizeof(unsigned char **));
  result.t_arr = t_arr;
  for (size_t i = 0; i < fund.len; i++) {
    t_arr[i] = calloc(cube.len - i - 1, sizeof(unsigned char *));
    for (size_t j = i + 1; j < cube.len; j++) {
      t_arr[i][j - i -1] = calloc((cube.len - j-  1) / 8 + 1, sizeof(unsigned char));
      if (init_value == MEM_LIST_TRUE) 
        for (size_t k = j+1; k < cube.len; k++) {
          tri_index indices = {i,j - i - 1,k - j - 1}; //Correct order as i < j < k
          SMI(t_arr,indices);
        }
    }
  }
  /*
   * Generate the sym_vertex_index. With this array we can convert
   * vertices to each of its 48 symmetries
  vert_index ** sym_vertex_index = malloc(cube.len * sizeof(vert_index *));
  for (unsigned short i = 0; i < cube.len; i++)
  {
    arr3 cur_pt;
    copyArr3(cur_pt, result.mem_fund.vert_from_index[i]);
    sym_vertex_index[i] = malloc(48 * sizeof(vert_index));
    for (int s = 0; s < 48; s++) {
      arr3 sym_pt;
      apply_symmetry(s,dim,cur_pt,sym_pt);
      sym_vertex_index[i][s] = vertex_to_index_fund(sym_pt, result.mem_fund.vert_to_index);
    }
  }
  result.mem_fund.sym_index = sym_vertex_index;

  /*
   * For each vertex we give the symmetry number needed to transform this
   * point into the fundamental domain
  int * vert_fund_sym = malloc(cube.len * sizeof(int));
  for (unsigned short i = 0; i < cube.len; i++) 
    for (int s = 0; s < 48; s++) 
      if (sym_vertex_index[i][s] < fund.len){ //This symmetry translates point i to the fund domain
        vert_fund_sym[i] = s;
        break;
      }

  result.mem_fund.vert_fund_sym = vert_fund_sym;

  free(cube.points);
  free(fund.points);
  return result;  
}
*/
