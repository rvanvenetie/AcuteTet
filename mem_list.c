#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include "vector.h"
#include "triangle.h"
#include "tetraeder.h"
#include "mem_list.h"
#include "combinations.h"

/*
 * Mem_list is a memory list used to store triangles inside the cube. 
 * 
 * Mem_list can have two options: 
 *  -a normal mem_list, in which it allocates stores for the first two
 *    dimensons and only allocates the third if an element in this dimension is needed.
 *  -a "fund" mem_list. Here only the triangles that reside in the fundamental area of 
 *    the cube are stored. The triangles are thus filtered by symmetry. To get 
 *    a symmetry version of the triangle, call triangle_to_index_fund. Indices
 *   not in the fundamental area are not initalized, thus one may call mem_list.t_arr[index] to check
 *   if index is inside the fundamental area.
 */
int index_list_from_file(tri_index_list * result, char * filename){
  FILE * stream;
  stream = fopen(filename, "rb");
  if (stream == NULL)
    return 0;
  if (fread(&result->dim, sizeof(arr3), 1, stream) < 1)
    return 0;
  if (fread(&result->len, sizeof(size_t), 1, stream) < 1)
    return 0;
  result->index_list = malloc(sizeof(tri_index) * result->len);
  if (fread(result->index_list, sizeof(tri_index), result->len, stream) < 1) {
    free(result->index_list);
    return 0;
  }
  fclose(stream);
  return 1;
}


int index_list_to_file(tri_index_list * list, char * filename) {
  FILE * stream;
  stream = fopen(filename, "wb");
  if (stream == NULL)
    return 0;
  //Write the struct to the file
  fwrite(list->dim, sizeof(arr3), 1, stream);
  fwrite(&list->len, sizeof(size_t), 1, stream);
  fwrite(list->index_list, sizeof(tri_index), list->len, stream);
  fclose(stream);
  return 1;  
}


//Writes mem_list to a file
int mem_list_to_file(tri_mem_list * list, char * filename) {
  FILE * stream;
  stream = fopen(filename, "wb");
  if (stream == NULL)
    return 0;
  //Write the struct to the file
  fwrite(list,sizeof(tri_mem_list), 1,stream);
  //Now the final structure is given by <c><last_dimension>
  //It consist of bytes. A 1 is followed by that dimensions value
  for (int i = 0; i < list->dim_size[0]; i++) {
    if (!list->t_arr[i])
      continue;
    for (int j = 0; j < list->dim_size[1]; j++) {
      if (list->t_arr[i][j]) { //Initalized
        fputc(1, stream); //Write 1 to indicate content is comming
        fwrite(list->t_arr[i][j],sizeof(unsigned char), list->dim_size[2], stream); //Write the contents
      } else { 
        fputc(0, stream); //No content is comming
      }
    }
  }
  fclose(stream);
  return 1;
}

//See write, but now in reverse
int mem_list_from_file(tri_mem_list * result, char * filename) {
  FILE * stream;
  stream = fopen(filename, "rb");
  if (stream == NULL)
    return 0;
  if (fread(result, sizeof(tri_mem_list), 1, stream) < 1)
    return 0;
  if (result->fund)
    *result = mem_list_init_fund(result->dim[0],MEM_LIST_FALSE);
  else
    *result = mem_list_init(result->dim,result->dim_size);
  for (int i = 0; i < result->dim_size[0]; i++) {
    if (!result->t_arr[i])
      continue;
    for (int j = 0; j < result->dim_size[1]; j++) {
      if (fgetc(stream)) {
        result->t_arr[i][j] = calloc(result->dim_size[2], sizeof(unsigned char));
        if (fread(result->t_arr[i][j], sizeof(unsigned char), result->dim_size[2], stream) < (unsigned int) result->dim_size[2])
          return 0;
      }
    }
  }
  fclose(stream);
  return 1;
}



tri_mem_list mem_list_init(arr3 dim, arr3 dim_size) { 
  tri_mem_list result = {NULL,NULL,
                        {dim_size[0],dim_size[1],dim_size[2]}, \
                        {(dim[1] + 1) * (dim[2] + 1),  (dim[2] + 1)}, \
                        {dim[0],dim[1],dim[2]},
                        0};
  unsigned char *** t_arr = calloc(dim_size[0], sizeof(unsigned char **));
 
  for (int i = 0; i < dim_size[0]; i++) 
    t_arr[i] = calloc(dim_size[1], sizeof(unsigned char *));  
  
  result.t_arr = t_arr;
  return result;
}


tri_mem_list mem_list_init_fund(int dim, int init_value) { 
  tri_mem_list result = {NULL,
                         NULL,
                        {(dim+1)*(dim+1)*(dim+1), (dim+1)*(dim+1)*(dim+1), (dim+1)*(dim+1)*(dim+1) / 8 + 1}, \
                        {(dim+1)*(dim+1),  (dim+1)}, \
                        {dim,dim,dim},
                        1};
  
  
    
  cube_points fund_points = gen_fund_points(dim);                      
  unsigned char *** t_arr = calloc(result.dim_size[0], sizeof(unsigned char **));
  result.t_arr = t_arr;
  
  /* Initalize first dimension */
  for (size_t i = 0; i < fund_points.len; i++) {
    int index = vertex_to_index(fund_points.points[i], result.dim_mult);
    t_arr[index] = calloc(result.dim_size[1], sizeof(unsigned char *));  
  }  
  /* Initalize second and third dimension */
  for (size_t i = 0; i < fund_points.len; i++) {
    int index = vertex_to_index(fund_points.points[i], result.dim_mult);
    for (int j = 0; j < result.dim_size[1]; j++) {
      t_arr[index][j] = calloc(result.dim_size[2], sizeof(unsigned char));
      if (init_value == MEM_LIST_TRUE) {
        if (index == j) //Same points are not valid triangles
          continue;
        for (int k = 0; k < result.dim_size[1]; k++) {
          if (index ==k || j == k)  //Double points -> not valid triangles
            continue;
          tri_index indices;
          vertices_unique_fund(index,j,k,&result,indices);
          if (indices[0] == index && indices[1] == j && indices[2] == k) {
            SMI(t_arr,indices);
          }
        }
      }
    }
  }
     
  
  /*
  size_t len = 0;
  Dindex * comb = combinations_list(result.dim_size[1],2,&len);
  for (size_t l = 0; l < len; l++) {//Loop through all combinations
    int j = comb[l*2];
    int k = comb[l*2+1];
    for (size_t i = 0; i < fund_points.len; i++){
      int index = vertex_to_index(fund_points.points[i],result.dim_mult);
      if (j == index || k == index)
        continue;
      tri_index indices;    
      vertices_unique_fund(index,j,k,&result,indices);
      
      if (!GMI(t_arr,indices))
        printf("%d,%d,%d\n", indices[0],indices[1],indices[2]);
      SMI(t_arr,indices);
    }
  }
  */
  unsigned short ** sym_vertex_index = malloc((result.dim_size[0]) * sizeof(unsigned short *));
  for (unsigned short i = 0; i < result.dim_size[0]; i++)
  {
    arr3 cur_pt;
    vertex_from_index(i, result.dim_mult, cur_pt);
    sym_vertex_index[i] = malloc(48 * sizeof(unsigned short));
    for (int s = 0; s < 48; s++) {
      arr3 sym_pt;
      apply_symmetry(s,dim,cur_pt,sym_pt);
      sym_vertex_index[i][s] = vertex_to_index(sym_pt, result.dim_mult);
    }
  }
  result.sym_index = sym_vertex_index;
  free(fund_points.points);
  return result;
}

void mem_list_free(tri_mem_list * list) {
  for (int i = 0; i < list->dim_size[0]; i++) { 
    if (!list->t_arr[i])
      continue;
    for (int j = 0; j < list->dim_size[1]; j++) 
      if (list->t_arr[i][j])
        free(list->t_arr[i][j]);
    free(list->t_arr[i]);
  }
  free(list->t_arr);
  
  for (unsigned short i = 0; i < list->dim_size[0]; i++)
    free(list->sym_index[i]);
  free(list->sym_index);
  
}
//Frees up rows that are entirely zero
void mem_list_clean(tri_mem_list * list) {
  for (int i = 0; i < list->dim_size[0]; i++) { 
    if (!list->t_arr[i])
      continue;
    for (int j = 0; j < list->dim_size[1]; j++) 
      if (list->t_arr[i][j]) {
        int row_zero = 1;
        for (int k = 0; k < list->dim_size[2]; k++) {
          if (list->t_arr[i][j][k]) {
            row_zero = 0;
            break;
          }
        }
        if (row_zero)  {
          free(list->t_arr[i][j]);
          list->t_arr[i][j] = NULL;
        }
      }
  
  }
}
/*
 * Gathers a list of all indices set to 1. Returns the amount of set indices.
 */
size_t mem_list_indices(tri_mem_list * list, tri_index_list * index_list) {
  size_t count = mem_list_count(list);
  index_list->dim[0] = list->dim[0];
  index_list->dim[1] = list->dim[1];
  index_list->dim[2] = list->dim[2];
  
  int c = 0;
  index_list->index_list = malloc(sizeof(tri_index) * count);
  
  for (int i = 0; i < list->dim_size[0]; i++) {
    if (!list->t_arr[i])
      continue;
    for (int j = 0; j < list->dim_size[1]; j++) 
      if (list->t_arr[i][j]) 
        for (int k = 0; k < list->dim_size[2] * 8; k++) {
          tri_index index = {i,j,k};
          if (GMI(list->t_arr, index)) { //Acute triangle with index i,j,k
            memcpy(index_list->index_list + c, &index, sizeof(tri_index));
            c++;
          }
        }         
  }
  index_list->len = c;
  return c;
}

triangle_list mem_list_to_triangle_list(tri_mem_list * list) {
  tri_index_list index_list;
  size_t count = mem_list_indices(list,&index_list);
  triangle_list result= {NULL, count, {list->dim[0],list->dim[1],list->dim[1]}}; 
  ptriangle t_arr = malloc(count * sizeof(triangle));
  for (size_t i = 0; i < count; i ++) {
    t_arr[i] = triangle_from_index(index_list.index_list[i], list->dim);
  }
  result.t_arr = t_arr;
  return result;
}
/*
 * Create a so-called mem_list from the triangle list.
 * Creates a sort of three-dimensional table.
 * First dimension is the "unique"-first vertex of the triangle.
 * Second dimension the second.
 * Third dimension is only created when needed. It consists of a bit array, enough space to hold 
 * the third vertex of unique representation.
 * 
 * TODO: Divide third dimension in blocks
 */
 
//Does symmetry shit!
tri_mem_list mem_list_from_triangle_list(triangle_list * list) {
  arr3 dim_size;
  dim_size[0] = (list->dim[0] + 1) * (list->dim[1] + 1) * (list->dim[2] + 1);
  dim_size[1] = dim_size[0];
  dim_size[2] = dim_size[1]/8 + 1;  
  tri_mem_list result = mem_list_init(list->dim,dim_size);
  //Now lets set all the nodes, exciting!
  for (int i = 0; i < list->len; i++) 
    mem_list_set_sym(&result, &list->t_arr[i]);
  return result;  
}

//Adds all the symmetries
tri_mem_list mem_list_from_index_list(tri_index_list * list) {
  arr3 dim_size;
  dim_size[0] = (list->dim[0] + 1) * (list->dim[1] + 1) * (list->dim[2] + 1);
  dim_size[1] = dim_size[0];
  dim_size[2] = dim_size[1]/8 + 1;  
  tri_mem_list result = mem_list_init(list->dim,dim_size);
  //Now lets set all the nodes, exciting!
  for (size_t i = 0; i < list->len; i++) {
    triangle cur_triang = triangle_from_index(list->index_list[i], list->dim);
    mem_list_set_sym(&result, &cur_triang);
  }
  return result;  
}

//Only adds triangles with one vertex in the fundamental domain (filters symmetries)
tri_mem_list mem_list_from_index_list_fund(tri_index_list * list) {
  tri_mem_list result = mem_list_init_fund(list->dim[0],MEM_LIST_FALSE);
  //Now lets set all the nodes, exciting!
  for (size_t i = 0; i < list->len; i++) {
    triangle cur_triang = triangle_from_index(list->index_list[i], list->dim);
    mem_list_set_sym_fund(&result, &cur_triang);
  }
  return result;    
  
}

size_t mem_list_memory(tri_mem_list * list) {
  size_t result = 0;
  result += list->dim_size[0] * sizeof(unsigned char **); //First dimension
//  result += list->dim_size[0] * list->dim_size[1] * sizeof(unsigned char *); //Second dimension
  
  for (int i = 0; i < list->dim_size[0]; i++) 
    if (list->t_arr[i]) {
      result += list->dim_size[1] * sizeof(unsigned char *); //Second dimension
      for (int j = 0; j < list->dim_size[1]; j++) 
        if (list->t_arr[i][j]) 
          result += (list->dim_size[2]) * sizeof(unsigned char);
    }
  return result;  
}

static const unsigned char BitsSetTable256[256] = 
{
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    B6(0), B6(1), B6(1), B6(2)
};

size_t mem_list_count(tri_mem_list * list) {
  size_t result = 0;
  for (int i = 0; i < list->dim_size[0]; i++) 
    if (list->t_arr[i])
      for (int j = 0; j < list->dim_size[1]; j++) 
        if (list->t_arr[i][j]) 
          for (int k = 0; k < list->dim_size[2]; k++) 
            result += BitsSetTable256[list->t_arr[i][j][k]];          
  return result;
}

void mem_list_clear(tri_mem_list * list, ptriangle triang){
  tri_index indices;
  triangle_to_index(triang, list->dim_mult, indices);
  CMI(list->t_arr, indices);
}

void mem_list_set(tri_mem_list * list, ptriangle triang)
{
  tri_index index;
  triangle_to_index(triang, list->dim_mult, index);
  if (!list->t_arr[index[0]])
    list->t_arr[index[0]] = calloc(list->dim_size[0], sizeof(unsigned char *));
  if (!EMI(list->t_arr,index)) //Not yet initalized
    list->t_arr[index[0]][index[1]] = calloc(list->dim_size[2], sizeof(unsigned char));
  //Set bit
  SMI(list->t_arr,index);
}  

int  mem_list_get(tri_mem_list * list, ptriangle triang) {
  tri_index index;
  triangle_to_index(triang, list->dim_mult, index);
  return GMI(list->t_arr, index);
}


void mem_list_clear_sym(tri_mem_list * list, ptriangle triang){
  for (int l =0; l < 48; l++) {
    triangle sym_triang;
    triangle_symmetry(triang, l, list->dim[0], &sym_triang);
    tri_index indices;
    triangle_to_index(&sym_triang, list->dim_mult, indices);
    CMI(list->t_arr, indices);
  }
}

void mem_list_set_sym(tri_mem_list * list, ptriangle triang) {
  for (int j=0; j < 48; j++) {
    triangle sym_triang;
    triangle_symmetry(triang, j, list->dim[0], &sym_triang);
    tri_index index;
    triangle_to_index(&sym_triang, list->dim_mult, index);
    if (!EMI(list->t_arr,index)) //Not yet initalized
      list->t_arr[index[0]][index[1]] = calloc(list->dim_size[2], sizeof(unsigned char));
    
    //Set bit
    SMI(list->t_arr,index);
  }  
}
void mem_list_clear_sym_fund(tri_mem_list * list, ptriangle triang){
  unsigned short * sym_v1, * sym_v2, * sym_v3;
  //Get pointers to the symmetry indexes of this triangle
  sym_v1 = list->sym_index[vertex_to_index(triang->vertices[0], list->dim_mult)];
  sym_v2 = list->sym_index[vertex_to_index(triang->vertices[1], list->dim_mult)];
  sym_v3 = list->sym_index[vertex_to_index(triang->vertices[2], list->dim_mult)];
  for (int s =0;s < 48; s++) {
    //triangle sym_triang;
    //triangle_symmetry(triang, l, list->dim[0], &sym_triang);
    //if (!triangle_to_index_fund(&sym_triang,list,index)) //This symmetry has no point inside the fund domain
    
    
    tri_index index;
    if (!vertices_unique_fund(sym_v1[s], sym_v2[s], sym_v3[s],list,index)) //No point in fund domain
      continue;    
    CMI(list->t_arr, index);
  }
}

//Only sets the symmetries of this triangle inside the fundamental domain
void mem_list_set_sym_fund(tri_mem_list * list, ptriangle triang) {
  unsigned short * sym_v1, * sym_v2, * sym_v3;
  //Get pointers to the symmetry indexes of this triangle
  sym_v1 = list->sym_index[vertex_to_index(triang->vertices[0], list->dim_mult)];
  sym_v2 = list->sym_index[vertex_to_index(triang->vertices[1], list->dim_mult)];
  sym_v3 = list->sym_index[vertex_to_index(triang->vertices[2], list->dim_mult)];
  for (int s=0; s < 48; s++) {
    //triangle sym_triang;
    //triangle_symmetry(triang, j, list->dim[0], &sym_triang);
    //if (!triangle_to_index_fund(&sym_triang,list,index)) //This symmetry has no point inside the fund domain
    
    tri_index index;
    if (!vertices_unique_fund(sym_v1[s], sym_v2[s], sym_v3[s],list,index)) //No point in fund domain
      continue;
    if (!EMI(list->t_arr,index)) //Not yet initalized
      list->t_arr[index[0]][index[1]] = calloc(list->dim_size[2], sizeof(unsigned char));
    //Set bit
    SMI(list->t_arr,index);
  }  
}
//Returns whether the given triangle is inside the mem_list. Translates it to a symmetry residing inside
//the fundamental domain and then checks. 
int mem_list_get_sym_fund(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3){
  unsigned short * sym_v1, * sym_v2, * sym_v3;
  //Get pointers to the symmetry indexes of this triangle
  sym_v1 = list->sym_index[vertex_to_index(v1, list->dim_mult)];
  sym_v2 = list->sym_index[vertex_to_index(v2, list->dim_mult)];
  sym_v3 = list->sym_index[vertex_to_index(v3, list->dim_mult)];
  for (int s=0; s < 48; s++) {
    //triangle sym_triang;
    //triangle_symmetry(&triang, j, list->dim[0], &sym_triang);
    //if (!vertices_to_index_fund(sym_triang.vertices[0], sym_triang.vertices[1], sym_triang.vertices[2],list,index)) //This symmetry has no point inside the fund domain
    
    tri_index index;
    if (!vertices_unique_fund(sym_v1[s], sym_v2[s], sym_v3[s],list,index)) //No point in fund domain
      continue;
    
    //Should be initalized!
    //if (!EMI(list->t_arr,index)) //Not yet initalized
      //return 0;
    return GMI(list->t_arr,index);
  } 
  return 0;
}
#define swap(x,y) {t=x;x=y;y=t;}
/*
 * Translates the vertices to indices for the three different dimensions of the
 * mem_list. Gets a unique representation for each combination of 3 vertices, by 
 * sorting the index list from low to high.
 */
void vertices_to_index(arr3 v1, arr3 v2, arr3 v3, arr2 dim_mult, tri_index indices){
  int t;
  indices[0] = vertex_to_index(v1, dim_mult);
  indices[1] = vertex_to_index(v2, dim_mult);
  indices[2] = vertex_to_index(v3, dim_mult);
  if (indices[0] > indices[1]) swap(indices[0],indices[1]);
  if (indices[1] > indices[2]) swap(indices[1],indices[2]);
  if (indices[0] > indices[1]) swap(indices[0],indices[1]);   
}

void triangle_to_index(ptriangle triang, arr2 dim_mult, tri_index indices) {
  int t;
  indices[0] = vertex_to_index(triang->vertices[0], dim_mult);
  indices[1] = vertex_to_index(triang->vertices[1], dim_mult);
  indices[2] = vertex_to_index(triang->vertices[2], dim_mult);
  if (indices[0] > indices[1]) swap(indices[0],indices[1]);
  if (indices[1] > indices[2]) swap(indices[1],indices[2]);
  if (indices[0] > indices[1]) swap(indices[0],indices[1]); 
}

/*
 * If none of the vertices belong to the fundamental area, the function returns 0.
 * Then it sorts the indices with the following rule: the first vertex is the lowest vertex INSIDE 
 * the fundamental area. The last two vertices are then sorted on low < high.
 * This gives a unique representation for triangles with one vertex in the fundamental area.
 */
int vertices_unique_fund(unsigned short idx1, unsigned short idx2, unsigned short idx3, tri_mem_list * mem_list, tri_index indices) {
  int t;
  indices[0] = idx1;
  indices[1] = idx2;
  indices[2] = idx3;
  if (!mem_list->t_arr[indices[0]] && !mem_list->t_arr[indices[1]] && !mem_list->t_arr[indices[2]])
    return 0;
  //First sort
  if (indices[0] > indices[1]) swap(indices[0],indices[1]);
  if (indices[1] > indices[2]) swap(indices[1],indices[2]);
  if (indices[0] > indices[1]) swap(indices[0],indices[1]);   
  
  //Set index[0] to be the lowest index that is also inside the fundamental domain
  if (mem_list->t_arr[indices[0]]) //Lowest index is also in the fundamental domain, current order is correct
    return 1;
  else if (mem_list->t_arr[indices[1]]) { //Second index, swap first and second
    swap(indices[0],indices[1]);
    return 1;
  }else if (mem_list->t_arr[indices[2]]) { //Last becomes first, second becomes third, first becomes second 
    swap(indices[0],indices[2]); //0 <-> 2, 
    swap(indices[2],indices[1]); //2 <-> 1
    return 1;
  } else
    return 0;  
}


int triangle_to_index_fund(ptriangle triang, tri_mem_list * mem_list, tri_index indices) {
  return vertices_to_index_fund(triang->vertices[0],triang->vertices[1], triang->vertices[2], mem_list,indices);
}
void vertex_from_index(unsigned short index, arr2 dim_mult,arr3 vertex) {
  int t;
  vertex[0] = index / dim_mult[0];
  t = index % dim_mult[0];
  vertex[1] = t / dim_mult[1];
  vertex[2] = t % dim_mult[1];    
}
triangle triangle_from_index(tri_index indices, arr3 dim) {
  triangle res;
  int YZ = (dim[1] + 1) * (dim[2] + 1);
  int Z  = (dim[2] + 1);
  int t;
  res.vertices[0][0] = indices[0] / YZ;
  t = indices[0] % YZ;
  res.vertices[0][1] = t / Z;
  res.vertices[0][2] = t % Z;
  res.vertices[1][0] = indices[1] / YZ;
  t = indices[1] % YZ;
  res.vertices[1][1] = t / Z;
  res.vertices[1][2] = t % Z;
  res.vertices[2][0] = indices[2] / YZ;
  t = indices[2] % YZ;
  res.vertices[2][1] = t / Z;
  res.vertices[2][2] = t % Z;
  return res;
}

int index_in_array(int index, int * indices, size_t len) {
  for (size_t i = 0; i < len; i++)
    if (index == indices[i])
      return 1;
  return 0;
}



void tetra_to_index(ptetra tet, arr2 dim_mult, tet_index indices) {
  int t;
  indices[0] = vertex_to_index(tet->vertices[0], dim_mult);
  indices[1] = vertex_to_index(tet->vertices[1], dim_mult);
  indices[2] = vertex_to_index(tet->vertices[2], dim_mult);
  indices[3] = vertex_to_index(tet->vertices[3], dim_mult);
  int swapped = 1;
  while (swapped) {
    swapped = 0;
    for (int i = 1; i < 4; i++)
      if (indices[i-1] > indices[i]) {
        swap(indices[i-1], indices[i])
        swapped = 1;
      }
  }
}
