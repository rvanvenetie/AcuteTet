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
 * Mem_list is a memory list used to store triangles inside the cube. A triangle
 * can be set in the list (TRUE) or can not be set in the list (FALSE). Each triangle thus has it's own bit. 
 * A triangle has three vertices. To store the triangles inside the mem_list
 * we convert these vertices (x,y,z) to an index (0,0,0) -> 0 and (dim, dim, dim) -> last_index.
 * These indices are then represented on a unique way, so that every triangle
 * has a unique representation in indices. 
 * The actual storage is done by creating a three-dimensional bit-array. Where the first
 * dimension is the first unique index of the triangle and the second dimension the second unique
 * index. Finally some bit-tricks are used to convert the third vertex to a bit in the third dimension.
 * 
 * Mem_list can be of three types:
 *  -a cube mem_list: Used if need to store every triangle in the cube.
 * 
 *  -a fund mem_list. Here only the triangles that have a vertex in the fundamental area of 
 *    the cube are stored. The triangles are thus filtered by symmetry. To check
 *    if a triangle is inside the mem_list, symmetries are applied until we have a representation of that
 *    triangle in the fundamental region. Indices
 *    not in the fundamental area are not initalized, thus one may call mem_list.t_arr[index] to check
 *    if index is inside the fundamental area.
 *  
 *  -a tet mem_list. Creates a mem_list that allows storage of triangles in the unit tetrahedron.
 * 
 * 
 * TODO: Convert system using a combinatorial number system.
 */

/*
 * Vertex to index conversion functions are given as macros in mem_list.h. Each function
 * converts a vertex (x,y,z) to an index for the mem_list.
 *  
 * Below are then given the unique representation functions. They represent the three indices in a unique way. 
 * This property is important as then we can represent each triangle in a unique way (vertices may have 6 different orders)
 */
 
#define swap(x,y) {t=x;x=y;y=t;}
/*
 * Returns the unique order for indices of the mem_list_cube: idx1 < idx2 < idx3.
 */ 
void indices_unique_cube(vert_index idx1, vert_index idx2, vert_index idx3, tri_index indices) {
  int t;
  indices[0] = idx1;
  indices[1] = idx2;
  indices[2] = idx3;
  if (indices[0] > indices[1]) swap(indices[0],indices[1]);
  if (indices[1] > indices[2]) swap(indices[1],indices[2]);
  if (indices[0] > indices[1]) swap(indices[0],indices[1]);
}

/*
 * Returns the unique order for indices of the mem_list_tet. First we sort the
 * indices so that idx1 < idx2 < idx3. Because of this order the unique index is then given by
 * idx1, idx2 - idx1, idx3 - idx2. See mem_list_init for more information of the storage class of type TET
 */
void indices_unique_tet(vert_index idx1, vert_index idx2, vert_index idx3, tri_index indices) {
  int t;
  indices[0] = idx1;
  indices[1] = idx2;
  indices[2] = idx3;
  if (indices[0] > indices[1]) swap(indices[0],indices[1]);
  if (indices[1] > indices[2]) swap(indices[1],indices[2]);
  if (indices[0] > indices[1]) swap(indices[0],indices[1]);
  indices[2] = indices[2] - indices[1] - 1;
  indices[1] = indices[1] - indices[0] - 1;  
}

/*
 * If none of the vertices belong to the fundamental area, the function returns 0.
 * To get a unique order it sorts the indices with the following rule: the first vertex is the lowest vertex INSIDE 
 * the fundamental area. The last two vertices are then sorted on low < high.
 * This gives a unique representation for triangles with one vertex in the fundamental area.
 */
int indices_unique_fund(vert_index idx1, vert_index idx2, vert_index idx3, tri_mem_list * mem_list, tri_index indices) {
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

/*
 * Vertices_to_index_* are macros defined in mem_list.h.
 * They convert three vertices to indices and then return the unique representation for the according
 * mem_list type.
 * 
 * Triangle_to_index_* are macros defined in mem_list.h. Wrappers for the above
 */


/*
 * vertex_from_index_*: Functions convert an index back to a vertex (x,y,z).
 * 
 * TODO: Implement _tet and _fund
 */
void vertex_from_index_cube(vert_index index, arr2 dim_mult,arr3 vertex) {
  int t;
  vertex[0] = index / dim_mult[0];
  t = index % dim_mult[0];
  vertex[1] = t / dim_mult[1];
  vertex[2] = t % dim_mult[1];    
}

triangle triangle_from_index_cube(tri_index indices, arr2 dim_mult) {
  triangle res;
  vertex_from_index_cube(indices[0], dim_mult, res.vertices[0]);
  vertex_from_index_cube(indices[1], dim_mult, res.vertices[1]);
  vertex_from_index_cube(indices[2], dim_mult, res.vertices[2]);
  return res;
}

/*
 * Remember that t_arr is a three dimensional array. Depening on the 
 * list type, the dimensions may not be fixed: the length of t_arr[0] may be
 * different than that of t_arr[1]. 
 * 
 * This function returns the size of the dimension dim:
 * - dim = 0 means length of t_arr
 * - dim = 1 means length of t_arr[idx1]
 * - dim = 2 means length of t_arr[idx1][idx2] (in BITS!!). 
 */  

int mem_list_dim_size(tri_mem_list * list, int dim, int idx1, int idx2) {
  switch (list->mode) {
    case MEM_LIST_CUBE:
      return list->mem_cube.dim_size[dim];
    case MEM_LIST_FUND:
      return list->mem_fund.dim_size;
    case MEM_LIST_TET:
      if (dim == 0)
        return list->mem_tet.tet_len; 
      else if (dim == 1)
        return list->mem_tet.tet_len - idx1 - 1;
      else  
        return list->mem_tet.tet_len - idx1 - idx2 - 2;
  }
  printf("WAT! NOT POSSIBLE");
  exit(0);
  return 0;
}



tri_mem_list mem_list_init_cube(arr3 dim) { 
  int dim_size = (dim[2] + 1) * (dim[1] + 1) * (dim[0] + 1);
  tri_mem_cube mem_cube = {{dim_size, dim_size, dim_size}, \
                           {(dim[1] + 1) * (dim[2] + 1),  (dim[2] + 1)}};
                           
  tri_mem_list result = {NULL,
                         {mem_cube},
                        {dim[0],dim[1],dim[2]},
                        MEM_LIST_CUBE};
                        
  unsigned char *** t_arr = calloc(dim_size, sizeof(unsigned char **));
 
  for (int i = 0; i < dim_size; i++) 
    t_arr[i] = calloc(dim_size, sizeof(unsigned char *));  
  
  result.t_arr = t_arr;
  return result;
}


tri_mem_list mem_list_init_tet(int dim, int init_value) { 
  tri_mem_list result;
  memset(&result,0, sizeof(tri_mem_list));
  result.dim[0] = dim;
  result.dim[1] = dim;
  result.dim[2] = dim;
  result.mode = MEM_LIST_TET;
  
  cube_points tet = gen_tet_points(dim, &result.mem_tet.vert_to_index);
  result.mem_tet.tet_len = tet.len;
  
  unsigned char *** t_arr = calloc(tet.len, sizeof(unsigned char **));
  result.t_arr = t_arr;
  
  for (size_t i = 0; i < tet.len; i++) {
    t_arr[i] = calloc(tet.len - i - 1, sizeof(unsigned char *));
    for (size_t j = i + 1; j < tet.len; j++) {
      t_arr[i][j - i -1] = calloc((tet.len - j-  1) / 8 + 1, sizeof(unsigned char));
      if (init_value == MEM_LIST_TRUE) 
        for (size_t k = j+1; k < tet.len; k++) {
          tri_index indices = {i,j - i - 1,k - j - 1}; //Correct order as i < j < k
          SMI(t_arr,indices);
        }
    }
  }
  free(tet.points);
  return result;
}



tri_mem_list mem_list_init_fund(int dim, int init_value) { 
  tri_mem_fund mem_fund = { (dim+1)*(dim+1)*(dim+1), \
                   {(dim+1)*(dim+1),  (dim+1)}, \
                   NULL};

  tri_mem_list result;
  memset(&result,0,sizeof(tri_mem_list));
  result.mode = MEM_LIST_FUND;
  result.mem_fund = mem_fund;
  result.dim[0] = dim; result.dim[1] = dim; result.dim[2] = dim;
  
  
    
  cube_points fund_points = gen_fund_points(dim);                      
  unsigned char *** t_arr = calloc(result.mem_fund.dim_size, sizeof(unsigned char **));
  result.t_arr = t_arr;
  
  /* Initalize first dimension */
  for (size_t i = 0; i < fund_points.len; i++) {
    int index = vertex_to_index_cube(fund_points.points[i], result.mem_fund.dim_mult);
    t_arr[index] = calloc(result.mem_fund.dim_size, sizeof(unsigned char *));  
  }  
  /* Initalize second and third dimension */
  for (size_t i = 0; i < fund_points.len; i++) {
    int index = vertex_to_index_cube(fund_points.points[i], result.mem_fund.dim_mult);
    for (int j = 0; j < result.mem_fund.dim_size; j++) {
      t_arr[index][j] = calloc(result.mem_fund.dim_size / 8 + 1, sizeof(unsigned char));
      if (init_value == MEM_LIST_TRUE) {
        if (index == j) //Same points are not valid triangles
          continue;
        for (int k = 0; k < result.mem_fund.dim_size; k++) {
          if (index ==k || j == k)  //Double points -> not valid triangles
            continue;
          tri_index indices;
          indices_unique_fund(index,j,k,&result,indices);
          if (indices[0] == index && indices[1] == j && indices[2] == k) {
            SMI(t_arr,indices);
          }
        }
      }
    }
  }
     
  unsigned short ** sym_vertex_index = malloc((result.mem_fund.dim_size) * sizeof(unsigned short *));
  for (unsigned short i = 0; i < result.mem_fund.dim_size; i++)
  {
    arr3 cur_pt;
    vertex_from_index_cube(i, result.mem_fund.dim_mult, cur_pt);
    sym_vertex_index[i] = malloc(48 * sizeof(unsigned short));
    for (int s = 0; s < 48; s++) {
      arr3 sym_pt;
      apply_symmetry(s,dim,cur_pt,sym_pt);
      sym_vertex_index[i][s] = vertex_to_index_cube(sym_pt, result.mem_fund.dim_mult);
    }
  }
  result.mem_fund.sym_index = sym_vertex_index;
  free(fund_points.points);
  return result;
}

void mem_list_free(tri_mem_list * list) {
  for (int i = 0; i < mem_list_dim_size(list,0,-1,-1); i++) {
    if (!list->t_arr[i])
      continue;
    for (int j = 0; j < mem_list_dim_size(list,1,i,-1); j++) 
      free(list->t_arr[i][j]);
    free(list->t_arr[i]);
  }
  free(list->t_arr);
  
  if (list->mode == MEM_LIST_FUND) {
    for (int i = 0; i < list->mem_fund.dim_size; i++)
      free(list->mem_fund.sym_index[i]);
    free(list->mem_fund.sym_index);
  }
}


/*
 * Returwns whether the row t_arr[i][j] is empty. Meaning that 
 * every triangle with unique index i,j is not in the mem_list
 */
int mem_list_row_empty(tri_mem_list * list,int i, int j) {
  for (int k = 0; k < (mem_list_dim_size(list,2,i,j)/ 8 + 1); k++)
    if (list->t_arr[i][j][k])
      return 0;
  return 1;
}

/*
 * Frees up rows that are entirely zero
 */
void mem_list_clean(tri_mem_list * list) {
  printf("mem_list_clean not implemented");
  /*
  for (int i = 0; i < list->dim_size[0]; i++) { 
    if (!list->t_arr[i])
      continue;
    for (int j = 0; j < list->dim_size[1]; j++) 
      if (list->t_arr[i][j] && mem_list_row_empty(list,i,j)) {
        free(list->t_arr[i][j]);
        list->t_arr[i][j] = NULL;
      }          
  }
  */
}




/*
 * Gathers a list of all the indices in MEM_LIST. Returns the total amount of triangles.
 */
size_t mem_list_indices(tri_mem_list * list, tri_index_list * index_list) {
  size_t count = mem_list_count(list);
  index_list->dim[0] = list->dim[0];
  index_list->dim[1] = list->dim[1];
  index_list->dim[2] = list->dim[2];
  
  int c = 0;
  index_list->index_list = malloc(sizeof(tri_index) * count);
  
  for (int i = 0; i < mem_list_dim_size(list,0,-1,-1); i++) {
    if (!list->t_arr[i])
      continue;
    for (int j = 0; j < mem_list_dim_size(list,1,i,-1); j++) 
      if (list->t_arr[i][j]) 
        for (int k = 0; k < mem_list_dim_size(list,2,i,j); k++) {
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

/*
 * Returns the triangles in MEM_LIST.
 */
 
triangle_list mem_list_to_triangle_list(tri_mem_list * list) {
  tri_index_list index_list;
  size_t count = mem_list_indices(list,&index_list);
  triangle_list result= {NULL, count, {list->dim[0],list->dim[1],list->dim[1]}}; 
  ptriangle t_arr = malloc(count * sizeof(triangle));
  
  if (!(list->mode == MEM_LIST_TET))
    for (size_t i = 0; i < count; i ++)
      t_arr[i] = triangle_from_index_cube(index_list.index_list[i], list->dim);
  
  result.t_arr = t_arr;
  return result;
}

static const unsigned char BitsSetTable256[256] = 
{
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    B6(0), B6(1), B6(1), B6(2)
};

/*
 * Returns the amount of triangles set in mem_list
 */
size_t mem_list_count(tri_mem_list * list) {
  size_t result = 0;
  
    for (int i = 0; i < mem_list_dim_size(list,0,-1,-1); i++)
      for (int j = 0; j < mem_list_dim_size(list,1,i,-1); j++)
        for (int k = 0; k < mem_list_dim_size(list,2,i,j) / 8 + 1; k++)
          result += BitsSetTable256[list->t_arr[i][j][k]];    
  return result;
}

/*
 * Returns the total amount of memory used by mem_list in BYTES
 */
size_t mem_list_memory(tri_mem_list * list) {
  size_t result = 0;
  
  result += mem_list_dim_size(list,0,-1,-1) * sizeof(unsigned char **); //First dimension
  for (int i = 0; i < mem_list_dim_size(list,0,-1,-1); i++) 
    if (list->t_arr[i]) {
      result += mem_list_dim_size(list,1,i,-1) * sizeof(unsigned char *); //Second dimension
      for (int j = 0; j < mem_list_dim_size(list,1,i,-1); j++) {
        if (list->t_arr[i][j]) 
          result += (mem_list_dim_size(list,2,i,j)/8 + 1 ) * sizeof(unsigned char);
      }
    }
  return result;  
}

/*
 * Stores the mem_list to filename.
 */
int mem_list_to_file(tri_mem_list * list, char * filename, int mode) {
  FILE * stream;
  stream = fopen(filename, "wb");
  if (stream == NULL)
    return 0;
  //Write the struct to the file
  if (fwrite(list,sizeof(tri_mem_list), 1,stream) < 1)
    return 0;


  //Now the final structure is given by <c><last_dimension>
  //It consist of bytes. A 1 is followed by that dimensions value
  for (int i = 0; i < mem_list_dim_size(list,0,-1,-1); i++) {
    if (!list->t_arr[i])
      continue;
    for (int j = 0; j < mem_list_dim_size(list,1,i,-1); j++) {
      if (list->t_arr[i][j] &&  //Initalized
         ((mode == MEM_LIST_SAVE_FULL) || //Save any way
         ((mode == MEM_LIST_SAVE_CLEAN) && !mem_list_row_empty(list,i,j)))) {
        fputc(1, stream); //Write 1 to indicate content is comming
        size_t write_count = mem_list_dim_size(list,2,i,j)/8 + 1;
        if (fwrite(list->t_arr[i][j],sizeof(unsigned char), write_count, stream) < write_count) //Write the contents
          return 0;
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
  
 
  for (int i = 0; i < mem_list_dim_size(result,0,-1,-1); i++) {
    if (!result->t_arr[i])
      continue;

    for (int j = 0; j < mem_list_dim_size(result,1,i,-1); j++) {
      size_t read_count = mem_list_dim_size(result,2,i,j);
      if (fgetc(stream)) {
        if (!result->t_arr[i][j])
          result->t_arr[i][j] = calloc(read_count, sizeof(unsigned char));
        if (fread(result->t_arr[i][j], sizeof(unsigned char), read_count, stream) < read_count)
          return 0;
      }
    }
  }
  fclose(stream);
  return 1;
}


/*
 * Remove a triangle from the mem_list of type CUBE
 */
void mem_list_clear_cube(tri_mem_list * list, ptriangle triang){
  tri_index indices;
  triangle_to_index_cube((*triang), list->mem_cube.dim_mult, indices);
  CMI(list->t_arr, indices);
}

/*
 * Add a triangle to the mem_list of type CUBE
 */
void mem_list_set_cube(tri_mem_list * list, ptriangle triang)
{
  tri_index index;
  triangle_to_index_cube((*triang), list->mem_cube.dim_mult, index);
  if (!list->t_arr[index[0]])
    list->t_arr[index[0]] = calloc(list->mem_cube.dim_size[0], sizeof(unsigned char *));
  if (!EMI(list->t_arr,index)) //Not yet initalized
    list->t_arr[index[0]][index[1]] = calloc(list->mem_cube.dim_size[2], sizeof(unsigned char));
  //Set bit
  SMI(list->t_arr,index);
}

/*
 * Return whether the given triangle is inside the mem_list of type CUBE
 */
int  mem_list_get_cube(tri_mem_list * list, ptriangle triang) {
  tri_index index;
  triangle_to_index_cube((*triang), list->mem_cube.dim_mult, index);
  return GMI(list->t_arr, index);
}

/*
 * Unset the triangle inside the mem_list of type CUBE + all of it's 48 symmetries
 */
void mem_list_clear_sym_cube(tri_mem_list * list, ptriangle triang){
  for (int l =0; l < 48; l++) {
    triangle sym_triang;
    triangle_symmetry(triang, l, list->dim[0], &sym_triang);
    tri_index indices;
    triangle_to_index_cube(sym_triang, list->mem_cube.dim_mult, indices);
    CMI(list->t_arr, indices);
  }
}

/*
 * Set the triangle inside the mem_list of type CUBE + all of it's 48 symmetries
 */
void mem_list_set_sym_cube(tri_mem_list * list, ptriangle triang) {
  for (int j=0; j < 48; j++) {
    triangle sym_triang;
    triangle_symmetry(triang, j, list->dim[0], &sym_triang);
    tri_index index;
    triangle_to_index_cube(sym_triang, list->mem_cube.dim_mult, index);
    if (!EMI(list->t_arr,index)) //Not yet initalized
      list->t_arr[index[0]][index[1]] = calloc(list->mem_cube.dim_size[2], sizeof(unsigned char));
    
    //Set bit
    SMI(list->t_arr,index);
  }  
}

/*
 * Unset the triangle inside the mem_list of type fund
 */
void mem_list_clear_fund(tri_mem_list * list, ptriangle triang){
  unsigned short * sym_v1, * sym_v2, * sym_v3;
  //Get pointers to the symmetry indexes of this triangle
  sym_v1 = list->mem_fund.sym_index[vertex_to_index_cube(triang->vertices[0], list->mem_fund.dim_mult)];
  sym_v2 = list->mem_fund.sym_index[vertex_to_index_cube(triang->vertices[1], list->mem_fund.dim_mult)];
  sym_v3 = list->mem_fund.sym_index[vertex_to_index_cube(triang->vertices[2], list->mem_fund.dim_mult)];
  for (int s =0;s < 48; s++) {
    //triangle sym_triang;
    //triangle_symmetry(triang, l, list->dim[0], &sym_triang);
    //if (!triangle_to_index_fund(&sym_triang,list,index)) //This symmetry has no point inside the fund domain
    
    
    tri_index index;
    if (!indices_unique_fund(sym_v1[s], sym_v2[s], sym_v3[s],list,index)) //No point in fund domain
      continue;    
    CMI(list->t_arr, index);
  }
}

/*
 * Set the triangle inside the mem_list of type fund
 */
void mem_list_set_fund(tri_mem_list * list, ptriangle triang) {
  unsigned short * sym_v1, * sym_v2, * sym_v3;
  //Get pointers to the symmetry indexes of this triangle
  sym_v1 = list->mem_fund.sym_index[vertex_to_index_cube(triang->vertices[0], list->mem_fund.dim_mult)];
  sym_v2 = list->mem_fund.sym_index[vertex_to_index_cube(triang->vertices[1], list->mem_fund.dim_mult)];
  sym_v3 = list->mem_fund.sym_index[vertex_to_index_cube(triang->vertices[2], list->mem_fund.dim_mult)];
  for (int s=0; s < 48; s++) {
    //triangle sym_triang;
    //triangle_symmetry(triang, j, list->dim[0], &sym_triang);
    //if (!triangle_to_index_fund(&sym_triang,list,index)) //This symmetry has no point inside the fund domain
    
    tri_index index;
    if (!indices_unique_fund(sym_v1[s], sym_v2[s], sym_v3[s],list,index)) //No point in fund domain
      continue;

    //Set bit
    SMI(list->t_arr,index);
  }  
}
/*
 * Returns whether the triangle with vertices v1,v2 and v3 is inside the mem_list of type FUND.
 */
int mem_list_get_fund(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3){
  unsigned short * sym_v1, * sym_v2, * sym_v3;
  //Get pointers to the symmetry indexes of this triangle
  sym_v1 = list->mem_fund.sym_index[vertex_to_index_cube(v1, list->mem_fund.dim_mult)];
  sym_v2 = list->mem_fund.sym_index[vertex_to_index_cube(v2, list->mem_fund.dim_mult)];
  sym_v3 = list->mem_fund.sym_index[vertex_to_index_cube(v3, list->mem_fund.dim_mult)];
  for (int s=0; s < 48; s++) {
    //triangle sym_triang;
    //triangle_symmetry(&triang, j, list->dim[0], &sym_triang);
    //if (!vertices_to_index_fund(sym_triang.vertices[0], sym_triang.vertices[1], sym_triang.vertices[2],list,index)) //This symmetry has no point inside the fund domain
    
    tri_index index;
    if (!indices_unique_fund(sym_v1[s], sym_v2[s], sym_v3[s],list,index)) //No point in fund domain
      continue;
    
    //Should be initalized!
    //if (!EMI(list->t_arr,index)) //Not yet initalized
      //return 0;
    return GMI(list->t_arr,index);
  } 
  return 0;
}

/*
 * Returns whether the triangle with vertices v1,v2 and v3 is inside the mem_list of type FUND.
 */
int mem_list_get_tet(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3){
  tri_index index;
  vertices_to_index_tet(v1,v2,v3, list->mem_tet.vert_to_index, index);
  return GMI(list->t_arr, index);
}


void mem_list_set_tet(tri_mem_list * list, ptriangle triang){
  tri_index index;
  vertices_to_index_tet(triang->vertices[0],triang->vertices[1],triang->vertices[2], list->mem_tet.vert_to_index, index);
  SMI(list->t_arr, index);
}

void mem_list_clear_tet(tri_mem_list * list, ptriangle triang){
  tri_index index;
  vertices_to_index_tet(triang->vertices[0],triang->vertices[1],triang->vertices[2], list->mem_tet.vert_to_index, index);
  CMI(list->t_arr, index);
}

/*
 * BELOW is not used code. Delete?
 */

void tetra_to_index_cube(ptetra tet, arr2 dim_mult, tet_index indices) {
  int t;
  indices[0] = vertex_to_index_cube(tet->vertices[0], dim_mult);
  indices[1] = vertex_to_index_cube(tet->vertices[1], dim_mult);
  indices[2] = vertex_to_index_cube(tet->vertices[2], dim_mult);
  indices[3] = vertex_to_index_cube(tet->vertices[3], dim_mult);
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


/*
 * Converts an list of indices to a mem_list of fund type.
 */
tri_mem_list mem_list_from_index_list_fund(tri_index_list * list) {
  tri_mem_list result = mem_list_init_fund(list->dim[0],MEM_LIST_FALSE);
  //Now lets set all the nodes, exciting!
  for (size_t i = 0; i < list->len; i++) {
    triangle cur_triang = triangle_from_index_cube(list->index_list[i], list->dim);
    mem_list_set_fund(&result, &cur_triang);
  }
  return result;
}


/*
 * Converts an list of indices to a mem_list of cube type.
 * Sets all the symmetries of the indices as well.
 */
tri_mem_list mem_list_from_index_list(tri_index_list * list) {
  tri_mem_list result = mem_list_init_cube(list->dim);
  //Now lets set all the nodes, exciting!
  for (size_t i = 0; i < list->len; i++) {
    triangle cur_triang = triangle_from_index_cube(list->index_list[i], list->dim);
    mem_list_set_sym_cube(&result, &cur_triang);
  }
  return result;  
}

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

