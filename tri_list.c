#include <string.h>
#include <stdio.h>
#include "triangle.h"
#include "tri_list.h"
/*
 * Assumes we have a sorted array. Uses binary search
 * to find the index of i in the array. Returns -1 if the
 * integer is not found
 */
int int_arr_idx(int_arr * arr, int i) {
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
  int dim_size = (dim +  1) * (dim + 1) * (dim + 1);
  result.dim = dim;
                        
  int_arr ** t_arr = malloc(sizeof(int_arr *) * dim_size);
 
  for (int i = 0; i < dim_size; i++)  {
    t_arr[i] = calloc((dim_size - i ) , sizeof(int_arr));
    if (init_value == MEM_LIST_TRUE)
      for (int j = 1; j < dim_size - i; j++){
        t_arr[i][j].len = (dim_size - j - i - 1);
        t_arr[i][j].p_arr = malloc(t_arr[i][j].len * sizeof(unsigned short));
        for (int k = 1; k < dim_size - j - i; k++)
          t_arr[i][j].p_arr[k-1] = k;
      }
  }
  result.t_arr = t_arr;
  return result;
}

void tri_list_insert(tri_list * list, ptriangle  triang) {
  tri_index idx;
  triangle_to_index_cube((*triang), list->dim, idx);
  int_arr * points = &list->t_arr[idx[0]][idx[1]];
  if (tri_list_contains(list, triang))
    return;

  unsigned short * new_arr = malloc((points->len + 1) * sizeof(unsigned short));
  int i = 0;
  while (i < points->len && points->p_arr[i] < idx[2] ) {
    new_arr[i] = points->p_arr[i];
    i++;
  }
 
  //Insert on place i
  new_arr[i] = idx[2];

  //Copy rest of the array
  memcpy(new_arr + i + 1, points->p_arr + i, (points->len - i) * sizeof(unsigned short));

  if (points->p_arr)
    free(points->p_arr);
  points->len++;
  points->p_arr = new_arr;
}

int tri_list_get(tri_list * list, arr3 v1, arr3 v2, arr3 v3) {
  tri_index idx;
  vertices_to_index_cube(v1,v2,v3,list->dim,idx);
  int_arr * points = &list->t_arr[idx[0]][idx[1]];
  return (int_arr_idx(points, idx[2]) != -1);
}
int tri_list_contains(tri_list * list, ptriangle  triang) {
  tri_index idx;
  triangle_to_index_cube((*triang), list->dim, idx);
  int_arr * points = &list->t_arr[idx[0]][idx[1]];
  return (int_arr_idx(points, idx[2]) != -1);
}

void tri_list_remove(tri_list * list, ptriangle  triang) {
  tri_index idx;
  triangle_to_index_cube((*triang), list->dim, idx);
  int_arr * points = &list->t_arr[idx[0]][idx[1]];
  int i = int_arr_idx(points,idx[2]);
  if (i < 0)
    return;
  unsigned short * new_arr = malloc((points->len - 1) * sizeof(unsigned short));
  //Copy everything before index i to the new_arr
  memcpy(new_arr, points->p_arr, i * sizeof(unsigned short));
  memcpy(new_arr + i, points->p_arr + i + 1, (points->len - i - 1) * sizeof(unsigned short));
  points->len--;
  free(points->p_arr);
  points->p_arr = new_arr;
}

void tri_list_free(tri_list * list) {
  int dim_size = (list->dim +  1) * (list->dim + 1) * (list->dim + 1);
 
  for (int i = 0; i < dim_size; i++) {
    for (int j = 0; j < dim_size - i; j++)
      free(list->t_arr[i][j].p_arr);
    free(list->t_arr[i]);
  }
  free(list->t_arr);
}
size_t tri_list_count(tri_list * list) {
  int dim_size = (list->dim +  1) * (list->dim + 1) * (list->dim + 1);
  int result = 0;
  for (int i = 0; i < dim_size; i++) {
    for (int j = 0; j < dim_size - i; j++) 
      result += list->t_arr[i][j].len;
  }
  return result;
}

int tri_list_to_file(tri_list * list, char * filename) {
  int dim_size = (list->dim +  1) * (list->dim + 1) * (list->dim + 1);
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
  int dim_size = (list->dim +  1) * (list->dim + 1) * (list->dim + 1);
  FILE * stream;
  stream = fopen(filename, "rb");
  if (stream == NULL)
    return 0;
    
  //Reads the entire struct from file
  if (fread(list, sizeof(tri_list), 1, stream) < 1)
    return 0;
  
  *list = tri_list_init(list->dim,MEM_LIST_FALSE);
 
  for (int i = 0; i < dim_size; i++) {
    for (int j = 0; j < dim_size - i; j++) {
      if (fread(&list->t_arr[i][j].len, sizeof(list->t_arr[i][j].len), 1, stream) < 1)
        return 0;
      list->t_arr[i][j].p_arr = malloc(list->t_arr[i][j].len  * sizeof(unsigned short));
      if (fread(list->t_arr[i][j].p_arr, sizeof(unsigned short), list->t_arr[i][j].len, stream) < (size_t) list->t_arr[i][j].len)
        return 0;
    }
  }
  fclose(stream);
  return 1;
}

