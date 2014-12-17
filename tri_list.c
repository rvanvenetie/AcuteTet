#include <string.h>
#include <stdio.h>
#include "triangle.h"
#include "tri_list.h"
#include "omp.h"

/*
 * Assumes we have a sorted array. Uses binary search
 * to find the index of i in the array. Returns -1 if the
 * integer is not found
 */
int int_arr_idx(p_int_arr arr, int i) {
  int lo = 0;
  int hi = arr->len - 1;

  while (lo <= hi) {
    int mid = (lo + hi) / 2;
    if (arr->p_arr[mid] == i)
      return mid;
    else if (arr->p_arr[mid] < i)
      lo = mid + 1;
    else
      hi = mid - 1;
  }
  return -1;
}

tri_list tri_list_init(int dim, int init_value) {
  tri_list result;
  int dim_size = tri_list_dim_size(dim);
  result.dim = dim;

  int_arr ** t_arr = malloc(sizeof(int_arr *) * dim_size);
 
  for (int i = 0; i < dim_size; i++)  {
    t_arr[i] = calloc((dim_size - i ) , sizeof(int_arr));
    if (init_value == MEM_LIST_TRUE) {
      for (int j = 1; j < dim_size - i; j++){
        t_arr[i][j].len      = (dim_size - j - i - 1);
	t_arr[i][j].data_len = (dim_size - j - i - 1);
        t_arr[i][j].p_arr = malloc(t_arr[i][j].len * sizeof(vert_index));
        for (int k = 1; k < dim_size - j - i; k++)
          t_arr[i][j].p_arr[k-1] = k;
      }
    }
  }
  result.t_arr = t_arr;
  return result;
}

void tri_list_insert(tri_list * list, ptriangle  triang, int resize) {
  if (tri_list_contains(list, triang)) //Already have this triangle..
    return;

  tri_index idx;
  triangle_to_index_cube((*triang), list->dim, idx);
  p_int_arr  points = &list->t_arr[idx[0]][idx[1]];

  /*
   * Two cases:
   * if we have enough memory, just find location and memmove the rest of the data
   * if we do not have enough memory, allocate a new block and copy data for the first part
   *   while looking for the index.
   */
  unsigned short i = 0; //Index of the data
  unsigned short * new_arr; //(New) location of the data
  if (points->data_len > points->len) // We have enough data
  {
    new_arr = points->p_arr; //We have enough room to store the new array in the old array
    //Find the location where to insert
    while (i < points->len && points->p_arr[i] < idx[2] ) 
      i++;

  } else {

    new_arr = malloc((points->len + 1) * sizeof(vert_index)); //New room needed
    //Find location and while at it, copy the first part of the array
    while (i < points->len && points->p_arr[i] < idx[2] ) {
      new_arr[i] = points->p_arr[i];
      i++;
    }
  }
  //Copy the data from [i .. end] to [i + 1 .. end+1]
  memmove(new_arr + i + 1, points->p_arr + i, (points->len - i) * sizeof(vert_index));
  //Insert on place i
  new_arr[i] = idx[2];
  //Edit length
  points->len++;


  if (new_arr != points->p_arr) { //New data array
    free(points->p_arr);
    points->p_arr = new_arr;
    points->data_len = points->len;
  }


  if (resize) {
    points->p_arr = realloc(points->p_arr, points->len * sizeof(vert_index));
    points->data_len = points->len;
  }
}

int tri_list_get(tri_list * list, arr3 v1, arr3 v2, arr3 v3) {
  tri_index idx;
  vertices_to_index_cube(v1,v2,v3,list->dim,idx);
  p_int_arr points = &list->t_arr[idx[0]][idx[1]];
  return (int_arr_idx(points, idx[2]) != -1);
}
int tri_list_contains(tri_list * list, ptriangle  triang) {
  tri_index idx;
  triangle_to_index_cube((*triang), list->dim, idx);
  p_int_arr  points = &list->t_arr[idx[0]][idx[1]];
  return (int_arr_idx(points, idx[2]) != -1);
}

void tri_list_remove(tri_list * list, ptriangle  triang, int resize) {
  tri_index idx;
  triangle_to_index_cube((*triang), list->dim, idx);
  p_int_arr points = &list->t_arr[idx[0]][idx[1]];
  int i = int_arr_idx(points,idx[2]);
  if (i < 0)
    return;
  //i holds the index of the point we want to remove. Copy data [i + 1 .. end] to [i .. end - 1]
  memmove(points->p_arr + i, points->p_arr + i + 1, (points->len - i -1) * sizeof(vert_index));
  //Decrease len
  points->len--;
  if (resize) 
  {
    points->p_arr = realloc(points->p_arr, points->len * sizeof(vert_index));
    points->data_len = points->len;
  }
  
}

void tri_list_empty(tri_list * list) {
  int dim_size = tri_list_dim_size(list->dim);
  for (int i = 0; i < dim_size; i++) 
    for (int j = 0; j < dim_size - i; j++) 
      list->t_arr[i][j].len = 0;
}

void tri_list_resize(tri_list * list) {
  int dim_size = tri_list_dim_size(list->dim);
  for (int i = 0; i < dim_size; i++) 
    for (int j = 0; j < dim_size - i; j++)  {
      list->t_arr[i][j].p_arr = realloc(list->t_arr[i][j].p_arr, list->t_arr[i][j].len * sizeof(vert_index));
      list->t_arr[i][j].data_len = list->t_arr[i][j].len;
    }
}
void tri_list_free(tri_list * list) {
  int dim_size = (list->dim +  1) * (list->dim + 1) * (list->dim + 1);
 
  for (int i = 0; i < dim_size; i++) {
    for (int j = 0; j < dim_size - i; j++){
      free(list->t_arr[i][j].p_arr);
    }
    free(list->t_arr[i]);
  }
  free(list->t_arr);
}

size_t tri_list_count(tri_list * list) {
  size_t dim_size = (list->dim +  1) * (list->dim + 1) * (list->dim + 1);
  size_t result = 0;
  for (size_t i = 0; i < dim_size; i++) {
    for (size_t j = 0; j < dim_size - i; j++) 
      result += list->t_arr[i][j].len;
  }
  return result;
}

size_t tri_list_memory(tri_list * list) {
  size_t dim_size = tri_list_dim_size(list->dim);
  size_t result = 0;

  result += sizeof(int_arr *) * dim_size;
  
 
  for (size_t i = 0; i < dim_size; i++)  {
    result += sizeof(int_arr) * (dim_size - i);
    for (size_t j = 1; j < dim_size - i; j++){
      result += sizeof(vert_index) * list->t_arr[i][j].data_len;
    }
  }
  return result;
}

int tri_list_to_file(tri_list * list, char * filename) {
  size_t dim_size = tri_list_dim_size(list->dim);
  FILE * stream;
  stream = fopen(filename, "wb");
  if (stream == NULL)
    return 0;
  //Write the struct to the file
  if (fwrite(list,sizeof(tri_list), 1,stream) < 1)
    return 0;

  for (int i = 0; i < dim_size; i++) {
    for (int j = 0; j < dim_size - i; j++) {
      fwrite(&list->t_arr[i][j].len, sizeof(list->t_arr[i][j].len), 1, stream);
      if (fwrite(list->t_arr[i][j].p_arr, sizeof(unsigned short), list->t_arr[i][j].len, stream) < (size_t) list->t_arr[i][j].len)
        return 0;
    }
  }
  fclose(stream);
  return 1;
}

int tri_list_from_file(tri_list * list, char * filename) {
  FILE * stream;
  stream = fopen(filename, "rb");
  if (stream == NULL)
    return 0;
    
  //Reads the entire struct from file
  if (fread(list, sizeof(tri_list), 1, stream) < 1)
    return 0;
  
  size_t dim_size = tri_list_dim_size(list->dim);
  *list = tri_list_init(list->dim,MEM_LIST_FALSE);
 
  for (size_t i = 0; i < dim_size; i++) {
    for (size_t j = 0; j < dim_size - i; j++) {
      if (fread(&list->t_arr[i][j].len, sizeof(list->t_arr[i][j].len), 1, stream) < 1)
        return 0;
      list->t_arr[i][j].data_len = list->t_arr[i][j].len;
      list->t_arr[i][j].p_arr = malloc(list->t_arr[i][j].len  * sizeof(unsigned short));
      if (fread(list->t_arr[i][j].p_arr, sizeof(unsigned short), list->t_arr[i][j].len, stream) < (size_t) list->t_arr[i][j].len)
        return 0;
    }
  }
  fclose(stream);
  return 1;
}
tri_list_old tri_list_init_old(int dim, int init_value) {
  tri_list_old result;
  int dim_size = (dim +  1) * (dim + 1) * (dim + 1);
  result.dim = dim;

  int_arr_old ** t_arr = malloc(sizeof(int_arr *) * dim_size);
 
  for (int i = 0; i < dim_size; i++)  {
    t_arr[i] = calloc((dim_size - i ) , sizeof(int_arr));
    if (init_value == MEM_LIST_TRUE) {
      for (int j = 1; j < dim_size - i; j++){
        t_arr[i][j].len      = (dim_size - j - i - 1);
        t_arr[i][j].p_arr = malloc(t_arr[i][j].len * sizeof(vert_index));
        for (int k = 1; k < dim_size - j - i; k++)
          t_arr[i][j].p_arr[k-1] = k;
      }
    }
  }
  result.t_arr = t_arr;
  return result;
}
int tri_list_old_from_file_to_new(tri_list * list, char * filename) {
  FILE * stream;
  stream = fopen(filename, "rb");
  if (stream == NULL)
    return 0;
    
  //Reads the entire struct from file
  if (fread(list, sizeof(tri_list), 1, stream) < 1)
    return 0;
  
  size_t dim_size = (list->dim +  1) * (list->dim + 1) * (list->dim + 1);
  *list = tri_list_init(list->dim,MEM_LIST_FALSE);
 
  for (size_t i = 0; i < dim_size; i++) {
    for (size_t j = 0; j < dim_size - i; j++) {
      size_t old_size;
      if (fread(&old_size, sizeof(old_size), 1, stream) < 1)
        return 0;
      list->t_arr[i][j].data_len = (unsigned short) old_size;
      list->t_arr[i][j].len      = (unsigned short) old_size;
      list->t_arr[i][j].p_arr = malloc(list->t_arr[i][j].len  * sizeof(unsigned short));
      if (fread(list->t_arr[i][j].p_arr, sizeof(unsigned short), list->t_arr[i][j].len, stream) < (size_t) list->t_arr[i][j].len)
        return 0;
    }
  }
  fclose(stream);
  return 1;
}

tri_list mem_list_to_tri_list(tri_mem_list * list) {
  if (list->mode != MEM_LIST_FUND) 
    return *(tri_list * )NULL;

  size_t dim_size = tri_list_dim_size(list->dim);
  tri_list result = tri_list_init(list->dim, MEM_LIST_FALSE);
  size_t i,j, cntr;
  unsigned short k;
  int sym_num; 
  tri_index cur_index, sym_index;
  unsigned short * tmp_array;
  arr3 v1,v2,v3;
  #pragma omp parallel private(j,k,i,v1,v2,v3, cntr, tmp_array, cur_index, sym_num, sym_index)
  {
    tmp_array = malloc(dim_size * sizeof(unsigned short));
    #pragma omp for schedule(dynamic, 20)
    for (i = 0; i < dim_size; i++) {

      vertex_from_index_cube(i,list->dim, v1); //Get point in cube
      cur_index[0] = vertex_to_index_fund(v1,list->mem_fund.vert_to_index); //Convert to index in mem_list
      sym_num = list->mem_fund.vert_fund_sym[cur_index[0]]; //Calculate symmetry needed
      cur_index[0] = list->mem_fund.sym_index[cur_index[0]][sym_num]; //Store symmetried point

      for ( j = 1; j < dim_size - i; j++) {
        cntr = 0;

        vertex_from_index_cube(i + j,list->dim, v2);
        cur_index[1] = vertex_to_index_fund(v2,list->mem_fund.vert_to_index);
        cur_index[1] = list->mem_fund.sym_index[cur_index[1]][sym_num];

        for (k = 1; k < dim_size - j - i; k++)
        {
          vertex_from_index_cube(i + j + k,list->dim, v3);
          cur_index[2] = vertex_to_index_fund(v3,list->mem_fund.vert_to_index);
          cur_index[2] = list->mem_fund.sym_index[cur_index[2]][sym_num];

          //Apply symmetry_number and return the unique_index of this transformed triangle.
          indices_unique_fund(cur_index[0],cur_index[1],cur_index[2],sym_index);

          //Sym_index now holds the index of the triangle given by v1,v2,v3 in FUND_MEM_LIST
          if (GMI(list->t_arr,sym_index))
            tmp_array[cntr++] = k;
        }
        result.t_arr[i][j].len = cntr;
        result.t_arr[i][j].data_len = cntr;
        result.t_arr[i][j].p_arr = malloc(cntr * sizeof(unsigned short));
        memcpy(result.t_arr[i][j].p_arr, tmp_array, cntr * sizeof(unsigned short));
      }
    }
  }
  return result;
}


/*
 *
 * Should be in a different file!
 */


data_list data_list_init(int dim, int mode, int init_value) {
  data_list result;
  result.mode = mode;
  switch (mode){
    case DATA_TRI_LIST:
      result.list = tri_list_init(dim,init_value);
      break;
    case DATA_MEM_LIST_FUND:
      result.mem_list = mem_list_init_fund(dim, init_value);
      break;
    case DATA_MEM_LIST_TET:
      result.mem_list = mem_list_init_tet(dim, init_value);
      break;
    case DATA_MEM_LIST_CUBE:
      result.mem_list = mem_list_init_cube(dim, init_value);
      break;
  }
  return result;
}
size_t data_list_count(data_list * list) {
  if (list->mode == DATA_TRI_LIST)
    return tri_list_count(&list->list);
  else
    return mem_list_count(&list->mem_list);
}
size_t data_list_memory(data_list * list) {
  if (list->mode == DATA_TRI_LIST)
    return tri_list_memory(&list->list);
  else
    return mem_list_memory(&list->mem_list);
}
void data_list_free(data_list * list) {
  if (list->mode == DATA_TRI_LIST)
    tri_list_free(&list->list);
  else
    mem_list_free(&list->mem_list);
}
int data_list_from_file(data_list * result, int mode, char * filename) {
  result->mode = mode;
  if (mode == DATA_TRI_LIST)
    return tri_list_from_file(&result->list, filename);
  else
    return mem_list_from_file(&result->mem_list, filename);

}
int data_list_to_file(data_list * list, char * filename, int mode) {
  if (mode == DATA_TRI_LIST)
    return tri_list_to_file(&list->list, filename);
  else
    return mem_list_to_file(&list->mem_list, filename, mode);
}
