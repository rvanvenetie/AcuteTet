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
  result.dim = dim;
  int dim_size = tri_list_dim_size(&result, 0, -1,-1);

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

int tri_list_insert(tri_list * list, ptriangle  triang, int resize) {
  if (tri_list_contains(list, triang)) //Already have this triangle..
    return 0;

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
  return 1;
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
  int dim_size = tri_list_dim_size(list, 0, -1,-1);
  for (int i = 0; i < dim_size; i++) 
    for (int j = 0; j < dim_size - i; j++) 
      list->t_arr[i][j].len = 0;
}

void tri_list_resize(tri_list * list) {
  int dim_size = tri_list_dim_size(list, 0, -1,-1);
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

void tri_list_validate(tri_list * list) {
  size_t dim_size = (list->dim +  1) * (list->dim + 1) * (list->dim + 1);
  for (size_t i = 0; i < dim_size; i++) 
    for (size_t j = 0; j < dim_size - i; j++) 
      for (size_t l = 1; l < list->t_arr[i][j].len; l++)
        if (list->t_arr[i][j].p_arr[l-1] >= list->t_arr[i][j].p_arr[l])
          printf("List is not valid!\n");
}
size_t tri_list_count(tri_list * list) {
  size_t dim_size = (list->dim +  1) * (list->dim + 1) * (list->dim + 1);
  size_t result = 0;
  for (size_t i = 0; i < dim_size; i++) 
    for (size_t j = 0; j < dim_size - i; j++) 
      if (list->t_arr[i][j].len) 
        result += list->t_arr[i][j].len;
  return result;
}

void tri_list_copy(tri_list * dest, tri_list * source) {
  if (dest->dim != source->dim) {
    tri_list_free(dest);
    *dest = tri_list_init(source->dim, MEM_LIST_FALSE);
  }

  int dim_size = tri_list_dim_size(source, 0, -1, -1);
  for (int i = 0; i < dim_size; i++)
    for (int j = 0; j < dim_size - i; j++) {
      if (source->t_arr[i][j].len > dest->t_arr[i][j].data_len) {
        free(dest->t_arr[i][j].p_arr);
        dest->t_arr[i][j].p_arr = malloc(source->t_arr[i][j].len * sizeof(vert_index));
        dest->t_arr[i][j].data_len = source->t_arr[i][j].len;
      }
      dest->t_arr[i][j].len = source->t_arr[i][j].len;
      memcpy(dest->t_arr[i][j].p_arr, source->t_arr[i][j].p_arr, dest->t_arr[i][j].len * sizeof(vert_index));
    }
}
size_t tri_list_memory(tri_list * list) {
  int dim_size = tri_list_dim_size(list, 0, -1,-1);
  size_t result = 0;

  result += sizeof(int_arr *) * dim_size;


  for (int i = 0; i < dim_size; i++)  {
    result += sizeof(int_arr) * (dim_size - i);
    for (int j = 1; j < dim_size - i; j++){
      result += sizeof(vert_index) * list->t_arr[i][j].data_len;
    }
  }
  return result;
}

int tri_list_dim_size(tri_list* list, int axis, int idx1, int idx2) {
  int dim_size = ((list->dim + 1) * (list->dim + 1) * (list->dim + 1));
  if (axis == 0)
    return dim_size;
  else if (axis == 1)
    return dim_size - idx1;
  else if (axis == 2) //If we have a sparse matrix, thirds axis might be empty
    return list->t_arr[idx1][idx2].len;
  return -1;
}
int tri_list_to_file(tri_list * list, char * filename) {
  int dim_size = tri_list_dim_size(list, 0, -1, -1);
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

  int dim_size = tri_list_dim_size(list,0,-1,-1);
  *list = tri_list_init(list->dim,MEM_LIST_FALSE);

  for (int i = 0; i < dim_size; i++) {
    for (int j = 0; j < dim_size - i; j++) {
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
int tri_list_update_from_file(tri_list * list, char * filename) {

  FILE * stream;
  stream = fopen(filename, "rb");
  if (stream == NULL)
    return 0;

  //Reads the entire struct from file
  if (fseek(stream, sizeof(tri_list), SEEK_CUR))
    return 0;

  int dim_size = tri_list_dim_size(list,0,-1,-1);

  for (int i = 0; i < dim_size; i++) {
    for (int j = 0; j < dim_size - i; j++) {
      //Read len parameter
      if (fseek(stream, sizeof(list->t_arr[i][j].len), SEEK_CUR))
        return 0;
      if (list->t_arr[i][j].data_len != list->t_arr[i][j].len) //We must read from file!
      {
        if (fread(list->t_arr[i][j].p_arr, sizeof(unsigned short), list->t_arr[i][j].len, stream) < (size_t) list->t_arr[i][j].len)
          return 0;
        list->t_arr[i][j].len = list->t_arr[i][j].data_len;
      } else
        if (fseek(stream, sizeof(unsigned short) * list->t_arr[i][j].len, SEEK_CUR))
          return 0;
    }
  }
  fclose(stream);
  return 1;

}
tri_list mem_list_to_tri_list(tri_mem_list * list) {
  if (list->mode != MEM_LIST_FUND) 
    return *(tri_list * )NULL;
  int dim_size = list->mem_fund.cube_len;
  tri_list result = tri_list_init(list->dim, MEM_LIST_FALSE);
  fprintf(stderr,"Size of tri_mem_list: %zu\n", mem_list_memory(list));
  fprintf(stderr,"Size of tri_list:     %zu\n", tri_list_memory(&result));
  int i,j;
  int s;
  unsigned short k;
  tri_index cur_idx, sym_idx;
  triangle cur_tri, sym_triang;
  arr3 * vert_from_index = list->mem_fund.vert_from_index;         
  omp_lock_t ** locks = malloc(sizeof(omp_lock_t *) * list->mem_fund.cube_len);
  //Initalize the locks
  for (i = 0; i < dim_size; i++){
    locks[i] = malloc(sizeof(omp_lock_t) * (dim_size));
    for (j = 0; j < dim_size - i; j++)
      omp_init_lock(&locks[i][j]);
  }

#pragma omp parallel for schedule(dynamic,5) shared(locks)\
  private(cur_idx, sym_idx, cur_tri, sym_triang,i,j,k,s)
  for (i = 0; i < mem_list_dim_size(list,0,-1,-1); i++)  {
    for (j = 0; j < mem_list_dim_size(list,1,i,-1); j++) {
      for (k = 0; k < mem_list_dim_size(list,2,i,j); k++) {
        cur_idx[0] = i;
        cur_idx[1] = j;
        cur_idx[2] = k;
        if (!GMI(list->t_arr, cur_idx))
          continue;
        cur_tri = triangle_from_index_fund(cur_idx, vert_from_index);

        for (s = 0; s < 48; s++) {
          triangle_symmetry(&cur_tri,s,list->dim,&sym_triang);
          triangle_to_index_cube((sym_triang), list->dim, sym_idx);

          omp_set_lock(&locks[sym_idx[0]][sym_idx[1]]);
          tri_list_insert(&result, &sym_triang, TRI_LIST_NO_RESIZE);
          omp_unset_lock(&locks[sym_idx[0]][sym_idx[1]]);
        }
      }
      free(list->t_arr[i][j]);
    }
    free(list->t_arr[i]);
    fprintf(stderr,"Thread %d is finnished with %d/%zu\n", omp_get_thread_num(), i, list->mem_fund.fund_len);
    if (i % 20 == 0)
      fprintf(stderr,"Amount of memory: %zu\n", tri_list_memory(&result)); 
  }
  for (i = 0; i < dim_size; i++){
    for (j = 0; j < dim_size -i ; j++)
      omp_destroy_lock(&locks[i][j]);
    free(locks[i]);
  }
  return result;
  /*
     arr3 v1,v2,v3;
     size_t cnt;
     tri_index cur_index, sym_index;
     unsigned short * tmp_array;
     int sym_num; 
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

  fprintf(stderr,"Thread %d finished doing %d\n",omp_get_thread_num(),i);
  }
  }
  */
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
      result.mem_list = mem_list_fund_init(dim, init_value, MEM_LIST_FUND);
      break;
    case DATA_MEM_LIST_FUND_SPARSE:
      result.mem_list = mem_list_fund_init(dim, init_value, MEM_LIST_FUND_SPARSE);
      break;
    case DATA_MEM_LIST_TET:
      result.mem_list = mem_list_tet_init(dim, init_value);
      break;
    case DATA_MEM_LIST_CUBE:
      result.mem_list = mem_list_cube_init(dim, init_value, MEM_LIST_CUBE_SPARSE);
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
int data_list_dim(data_list * list) {
  if (list->mode == DATA_TRI_LIST)
    return list->list.dim;
  else
    return list->mem_list.dim;
}
int data_list_dim_size(data_list * list, int dim, int idx1, int idx2) {
  if (list->mode == DATA_TRI_LIST)
    return tri_list_dim_size(&list->list, dim, idx1, idx2);
  else
    return mem_list_dim_size(&list->mem_list, dim, idx1,idx2);
}

int data_list_contains(data_list * list, ptriangle  triang) {
  if (list->mode == DATA_TRI_LIST)
    return tri_list_contains(&list->list, triang);
  else
    return mem_list_cube_contains(&list->mem_list, triang);

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
