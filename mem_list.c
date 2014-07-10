#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include "vector.h"
#include "triangle.h"
#include "tetraeder.h"
#include "mem_list.h"

/*
 * Mem_list is a memory list used to store a set of facets inside the cube. A triangle
 * can be set in the list (TRUE) or can not be set in the list (FALSE). Each possible triangle thus has it's own bit. 
 * A triangle has three vertices: v1,v2 and v3. To store the triangles we convert these vertices (points in 3D, thus array of 3 integers)
 * to indices. The vertices of a triangle are thus represented by three indices.
 * These indices are then represented on a unique way, so that every triangle
 * has a unique representation in memory_list.
 *  
 * The actual storage is done by creating a three-dimensional bit-array. Where the first
 * axis is the first unique index of the triangle and the second axis the second unique
 * index. Finally some bit-tricks are used to convert the third index to a bit in the third dimension.
 * 
 * Mem_list can be of three types:
 *  -a cube mem_list: Used if need to store every (possible) triangle in the cube.
 * 
 *  -a fund mem_list: Here only the triangles that have a vertex in the fundamental area of 
 *    the cube are stored. One may apply symmetries on a `normal' triangle, to find it's fundamental
 *    representation.
 *  
 *  -a tet mem_list: Used to store every possible triangle in the unit tetrahedron.
 * 
 * See the SMI/GMI/CMI macro's in mem_list.h for the setting etc. of triangles.
 * 
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
 * Turns 3 indices into the unique number used in mem_list.
 * First sorts the indices so idx1 < idx2 < idx3. Because this is the unique order
 * array indices index is then given by (idx1, idx2 - idx1, idx3 - idx2).
 */ 
void indices_unique_cube(vert_index idx1, vert_index idx2, vert_index idx3, tri_index indices) {
  int t;
  indices[0] = idx1;
  indices[1] = idx2;
  indices[2] = idx3;
  if (indices[0] > indices[1]) swap(indices[0],indices[1]);
  if (indices[1] > indices[2]) swap(indices[1],indices[2]);
  if (indices[0] > indices[1]) swap(indices[0],indices[1]);
  indices[2] = indices[2] - indices[1];
  indices[1] = indices[1] - indices[0];  
}

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

void indices_unique_fund(vert_index idx1, vert_index idx2, vert_index idx3, tri_index indices) {
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
 * Vertices_to_index_* are macros defined in mem_list.h. They convert
 * vertices to index numbers, which then are used in the functions above.
 */

/*
 * Generates a 3D array that converts a point to it's index. ie vertex_index[x][y][z] = idx of vertex (x,y,z)
 * Returns total amount of points.
 */
void gen_vertex_to_index_tet(int dim, vert_index **** vertex_index, int *tet_len) {
  int c = 0;
  *vertex_index = malloc((dim+1) * sizeof(unsigned short **));
  for (int x = 0; x <= dim; x++) {
    (*vertex_index)[x] = malloc( (dim-x + 1) * sizeof(unsigned short *));
    for (int y = 0; y <= dim - x; y++) {
      (*vertex_index)[x][y] = malloc( (dim - x - y + 1) * sizeof(unsigned short));      
      for (int z = 0; z <= dim - x -y; z++) {
        (*vertex_index)[x][y][z] = c;
        c++;
      }
    }
  }
  *tet_len = c;  
}

/*
 * Generates a 3D array for converting a point to it's index. First stores all the points in the
 * fundamental domain (0..fund_len) then the rest of the points in the cube.
 */
void gen_vertex_to_index_fund(int dim, vert_index_array * vertex_index, size_t * fund_len, size_t * cube_len) {
  cube_points fund_pts = gen_fund_points(dim);
  cube_points cube_pts = gen_cube_points(dim);
  *fund_len = fund_pts.len;
  *cube_len = cube_pts.len;
  
  //Allocate the 3D array, initalize to USHRT_MAX
  *vertex_index = malloc((dim+1) * sizeof(unsigned short **));
  for (int x = 0; x <= dim; x++) {
    (*vertex_index)[x] = malloc((dim + 1) * sizeof(unsigned short *));
    for (int y = 0; y <= dim; y++) {
      (*vertex_index)[x][y] = malloc((dim + 1) * sizeof(unsigned short));
      memset((*vertex_index)[x][y], USHRT_MAX, (dim + 1) * sizeof(unsigned short));
    }
  }
  //Fill every entry with the number of this vertex
  int c = 0;
  //First do the fundamental domain
  for (size_t i = 0; i < fund_pts.len; i ++)
  {
    (*vertex_index)[fund_pts.points[i][0]][fund_pts.points[i][1]][fund_pts.points[i][2]] = c;
    c++;
  }
  //Now fill the entire cube - fundamental domain in the rest. 
  for (size_t i = 0; i < cube_pts.len; i++)
    if ((*vertex_index)[cube_pts.points[i][0]][cube_pts.points[i][1]][cube_pts.points[i][2]] == USHRT_MAX) //Point is not inside the fund domain
    {
      (*vertex_index)[cube_pts.points[i][0]][cube_pts.points[i][1]][cube_pts.points[i][2]] = c;
      c++;
    }
  free(cube_pts.points);
  free(fund_pts.points);
}

/*
 * Generates an array that does the revers from above, ie index to point.
 */
void gen_vertex_from_index_fund(int dim, arr3 ** index_vertex, vert_index_array vertex_index) {
  size_t len = (dim+1) * (dim + 1) * (dim + 1); //Amount of points
  (*index_vertex) = malloc(len * sizeof(arr3));
  for (int x = 0; x <= dim; x++)
    for (int y = 0; y <= dim; y++)
      for (int z = 0; z <= dim; z++) {
        vert_index index = vertex_index[x][y][z];
        (*index_vertex)[index][0] = x;
        (*index_vertex)[index][1] = y;
        (*index_vertex)[index][2] = z;
     }
}

/*
 * vertex_from_index_*: Functions convert an index back to a vertex (x,y,z).
 */
 
//Conversion for points in cube
void vertex_from_index_cube(vert_index index, int dim ,arr3 vertex) {
  int t;
  vertex[0] = index / ((dim+1)*(dim+1));
  t = index % ((dim+1)*(dim+1));
  vertex[1] = t / (dim+1);
  vertex[2] = t % (dim+1);    
}

triangle triangle_from_index_cube(tri_index indices, int dim) {
  triangle res;
  vertex_from_index_cube(indices[0], dim, res.vertices[0]);
  vertex_from_index_cube(indices[1], dim, res.vertices[1]);
  vertex_from_index_cube(indices[2], dim, res.vertices[2]);
  return res;
}
triangle triangle_from_index_fund(tri_index indices,arr3 * index_vertex) {
  triangle res;
  copyArr3(res.vertices[0], index_vertex[indices[0]]);
  copyArr3(res.vertices[1], index_vertex[indices[0] + indices[1] + 1]);
  copyArr3(res.vertices[2], index_vertex[indices[0] + indices[1] + indices[2] + 2]);
  return res;
}
/*
 * Remember that t_arr is a three dimensional array. Depening on the 
 * list type, the dimensions may not be fixed: the length of t_arr[0] may be
 * different than that of t_arr[1]. 
 * 
 * This function returns the size of the dimension dim:
 * - axis = 0 means length of t_arr
 * - axis = 1 means length of t_arr[idx1]
 * - axis = 2 means length of t_arr[idx1][idx2] (in BITS!!). 
 */  

int mem_list_dim_size(tri_mem_list * list, int axis, int idx1, int idx2) {
  switch (list->mode) {
    case MEM_LIST_CUBE:
      if (axis == 0)
        return list->mem_cube.dim_size;
      else if (axis == 1)
        return list->mem_cube.dim_size - idx1;
      else if (axis == 2)
        return list->mem_cube.dim_size - idx1 - idx2 ;
    case MEM_LIST_fund:
      if (axis == 0)
        return list->mem_fund.fund_len;
      else if (axis == 1)
        return list->mem_fund.cube_len - idx1 - 1;
      else if (axis == 2)
        return list->mem_fund.cube_len - idx1 - idx2 - 2;
            
    case MEM_LIST_TET:
      if (axis == 0)
        return list->mem_tet.tet_len; 
      else if (axis == 1)
        return list->mem_tet.tet_len - idx1 - 1;
      else  
        return list->mem_tet.tet_len - idx1 - idx2 - 2;
  }
  printf("WAT! NOT POSSIBLE");
  exit(0);
  return 0;
}

/*
 * Initalize a memory of type MODE for object with dimension dim.
 * Init_value is  MEM_LIST_FALSE or MEM_LIST_TRUE. Setting this
 * can initalize to TRUE makes it set all the triangles in the memory list.
 */
tri_mem_list mem_list_init(int dim, int mode, int init_value) {
  printf("Initalizing memory_list: %d, %d, %d\n", dim, mode, init_value);
  switch (mode) {
    case MEM_LIST_CUBE:
      return mem_list_init_cube(dim,init_value);
    case MEM_LIST_fund:
      return mem_list_init_fund(dim, init_value);
    case MEM_LIST_TET:
      return mem_list_init_tet(dim, init_value);
  }
  //return NULL;
}

/*
 * Initalize memory_list of type cube
 */
tri_mem_list mem_list_init_cube(int dim,int init_value) { 
  int dim_size = (dim + 1) * (dim + 1) * (dim + 1);
  tri_mem_cube mem_cube = {dim_size};
                           
  tri_mem_list result = {NULL,
                         {mem_cube},
                        dim,
                        MEM_LIST_CUBE};
                        
  unsigned char *** t_arr = calloc(dim_size, sizeof(unsigned char **));
 
  for (int i = 0; i < dim_size; i++) {
    t_arr[i] = calloc(dim_size - i , sizeof(unsigned char *));
    for (int j = i; j < dim_size; j++) {
      t_arr[i][j - i] = calloc((dim_size - j) / 8 + 1, sizeof(unsigned char));
      if (init_value == MEM_LIST_TRUE) 
        for (int k = j+1; k < dim_size; k++) {
          if (i < j && j < k) {
            tri_index tmp_index = {i,j - i ,k - j };
            SMI(t_arr,tmp_index);
          }
        }
    }
  }

  
  result.t_arr = t_arr;
  return result;
}


tri_mem_list mem_list_init_tet(int dim, int init_value) { 
  tri_mem_list result;
  memset(&result,0, sizeof(tri_mem_list));
  result.dim = dim;
  result.mode = MEM_LIST_TET;
  
  cube_points tet = gen_tet_points(dim);
  //Generate vertex to index array
  gen_vertex_to_index_tet(dim, &result.mem_tet.vert_to_index, &result.mem_tet.tet_len);
  
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
  tri_mem_list result;
  memset(&result,0,sizeof(tri_mem_list));
  result.mode = MEM_LIST_fund;
  result.dim = dim; 
    
  //Generate vertex to index array  
  gen_vertex_to_index_fund(dim, 
                        &result.mem_fund.vert_to_index,
                        &result.mem_fund.fund_len,
                        &result.mem_fund.cube_len);
  gen_vertex_from_index_fund(dim, &result.mem_fund.vert_from_index, result.mem_fund.vert_to_index);

  //Cache values, easier to read..
  size_t fund_len = result.mem_fund.fund_len;
  size_t cube_len = result.mem_fund.cube_len;
  unsigned char *** t_arr = calloc(fund_len, sizeof(unsigned char **));
  result.t_arr = t_arr;
  
  /* Initalize the 3D bit array.
   * First dimension only holds points in fundamental domain (the first fund_len points)
   * Second dimension then holds all the points bigger than first dimension (also outside fund domain)
   * Third dimension holds all points bigger than the third dimension
   */
  for (size_t i = 0; i < fund_len; i++) {
    t_arr[i] = calloc(cube_len - i - 1, sizeof(unsigned char *));
    for (size_t j = i + 1; j < cube_len; j++) {
      t_arr[i][j - i -1] = calloc((cube_len - j-  1) / 8 + 1, sizeof(unsigned char));
      if (init_value == MEM_LIST_TRUE) 
        for (size_t k = j+1; k < cube_len; k++) {
          tri_index indices = {i,j - i - 1,k - j - 1}; //Correct order as i < j < k
          SMI(t_arr,indices);
        }
    }
  }
  vert_index ** sym_vertex_index = malloc(cube_len * sizeof(vert_index *));
  for (unsigned short i = 0; i < cube_len; i++)
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
  
  int * vert_fund_sym = malloc(cube_len * sizeof(int));
  for (unsigned short i = 0; i < cube_len; i++) 
    for (int s = 0; s < 48; s++) 
      if (sym_vertex_index[i][s] < fund_len){ //This symmetry translates point i to the fund domain
        vert_fund_sym[i] = s;
        break;
      }
  
  result.mem_fund.vert_fund_sym = vert_fund_sym;
  
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
  
  if (list->mode == MEM_LIST_fund) {
    vert_index ** sym_index;
    sym_index = list->mem_fund.sym_index;
      
    for (int i = 0; i < mem_list_dim_size(list,0,-1,-1); i++)
      free(sym_index[i]);
    free(sym_index);
    
    free(list->mem_fund.vert_fund_sym);
  }
}


/*
 * returns whether the row t_arr[i][j] is empty. Meaning that 
 * every triangle with unique index i,j is not in the mem_list
 */
int mem_list_row_empty(tri_mem_list * list,int i, int j) {
  for (int k = 0; k < (mem_list_dim_size(list,2,i,j)/ 8 + 1); k++)
    if (list->t_arr[i][j][k])
      return 0;
  return 1;
}

/*
 * Gathers a list of all the indices in MEM_LIST. Returns the total amount of triangles.
 * (This thus returns a list of all the indices of the triangles set in the memory_list.
 */
size_t mem_list_indices(tri_mem_list * list, tri_index_list * index_list) {
  size_t count = mem_list_count(list);
  index_list->dim = list->dim;
  
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
 * Returns a list of triangles set in the memory_list
 */
 
triangle_list mem_list_to_triangle_list(tri_mem_list * list) {
  tri_index_list index_list;
  size_t count = mem_list_indices(list,&index_list);
  triangle_list result= {NULL, count, list->dim}; 
  ptriangle t_arr = malloc(count * sizeof(triangle));
  
  if (list->mode == MEM_LIST_fund) {
    for (size_t i = 0; i < count; i ++) {
      printf("%d,%d,%d\n", index_list.index_list[i][0],index_list.index_list[i][1],index_list.index_list[i][2]);
      t_arr[i] = triangle_from_index_fund(index_list.index_list[i], list->mem_fund.vert_from_index); 
      print_triangle(t_arr + i);
    }
  }
  result.t_arr = t_arr;
  return result;
}

//Source: https://graphics.stanford.edu/~seander/bithacks.html

static const unsigned char BitsSetTable256[256] = 
{
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    B6(0), B6(1), B6(1), B6(2)
};

/*
 * Returns the amount of triangles set in mem_list.
 */
size_t mem_list_count(tri_mem_list * list) {
  size_t result = 0;
  
    for (int i = 0; i < mem_list_dim_size(list,0,-1,-1); i++)
     if (list->t_arr[i])
        for (int j = 0; j < mem_list_dim_size(list,1,i,-1); j++)
          if (list->t_arr[i][j])
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
 * Stores the mem_list to file. Mode can be either:
 * - MEM_LIST_SAVE_FULL: Save the entire mem_list to the file
 * - MEM_LIST_SAVE_CLEAN: Only save the rows that are not entirely zero to the file
 */
int mem_list_to_file(tri_mem_list * list, char * filename, int mode) {
  FILE * stream;
  stream = fopen(filename, "wb");
  if (stream == NULL)
    return 0;
  //Write the struct to the file
  if (fwrite(list,sizeof(tri_mem_list), 1,stream) < 1)
    return 0;

  //Saves to file by looping over first two axis.
  //Then it saves <c><last_dimension>
  //Here c indicates whether the row is saved (0 can occur if MEM_LIST_SAVE_CLEAN and row is empty)
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

//Loads mem_list_from_file, the above function in reverse.
int mem_list_from_file(tri_mem_list * result, char * filename) {
  FILE * stream;
  stream = fopen(filename, "rb");
  if (stream == NULL)
    return 0;
    
  //Reads the entire struct from file
  if (fread(result, sizeof(tri_mem_list), 1, stream) < 1)
    return 0;
  
  //Initalize to FALSE, as we are only going to set the triangles from the file
  *result = mem_list_init(result->dim, result->mode, MEM_LIST_FALSE);
 
  for (int i = 0; i < mem_list_dim_size(result,0,-1,-1); i++) {
    if (!result->t_arr[i])
      continue;

    for (int j = 0; j < mem_list_dim_size(result,1,i,-1); j++) {
      size_t read_count = mem_list_dim_size(result,2,i,j) / 8 + 1;
      if (fgetc(stream)) { //if fgetc = 1
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
 * Here folow the actual setting/getting/clearing functions for the mem_list.
 * All depending on the type of memory_list
 */

/*
 * Remove a triangle from the mem_list of type CUBE
 */
void mem_list_clear_cube(tri_mem_list * list, ptriangle triang){
  tri_index indices;
  triangle_to_index_cube((*triang), list->dim, indices);
  CMI(list->t_arr, indices);
}

/*
 * Add a triangle to the mem_list of type CUBE
 */
void mem_list_set_cube(tri_mem_list * list, ptriangle triang)
{
  tri_index index;
  triangle_to_index_cube((*triang), list->dim, index);
  SMI(list->t_arr,index);
}

/*
 * Return whether the given triangle is inside the mem_list of type CUBE
 */
int  mem_list_get_cube(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3) {
  tri_index index;
  vertices_to_index_cube(v1,v2,v3,list->dim,index);
  return GMI(list->t_arr, index);
}


/*
 * Clear triangle and all of it's symmetries from the mem_list
 */
void mem_list_clear_cube_sym(tri_mem_list * list, ptriangle triang){
  for (int i = 0; i < 48; i++){
    triangle sym_triang;
    triangle_symmetry(triang,i,list->dim,&sym_triang);
    tri_index indices;
    triangle_to_index_cube(sym_triang, list->dim, indices);
    CMI(list->t_arr, indices);
  }
}


/*
 * Clear triangle from the mem_list of type FUND. Does this by applying all symmetries on the vertices
 * and clearing every of the 48 symmetries.
 */


void mem_list_clear_fund(tri_mem_list * list, ptriangle triang){
  unsigned short * sym_v1, * sym_v2, * sym_v3;
  //Get pointers to the symmetry indexes of this triangle
  sym_v1 = list->mem_fund.sym_index[vertex_to_index_fund(triang->vertices[0], list->mem_fund.vert_to_index)];
  sym_v2 = list->mem_fund.sym_index[vertex_to_index_fund(triang->vertices[1], list->mem_fund.vert_to_index)];
  sym_v3 = list->mem_fund.sym_index[vertex_to_index_fund(triang->vertices[2], list->mem_fund.vert_to_index)];
  for (int s=0; s < 48; s++) {
    tri_index index;
    if (triangle_in_fund(sym_v1[s], sym_v2[s], sym_v3[s], list->mem_fund.fund_len))
    {
      indices_unique_fund(sym_v1[s],sym_v2[s],sym_v3[s], index);
      CMI(list->t_arr, index);
    }

  } 
}

/*
 * Does the same as the above function, but now sets a triangle
 */
void mem_list_set_fund(tri_mem_list * list, ptriangle triang) {
  unsigned short * sym_v1, * sym_v2, * sym_v3;
  //Get pointers to the symmetry indexes of this triangle
  sym_v1 = list->mem_fund.sym_index[vertex_to_index_fund(triang->vertices[0], list->mem_fund.vert_to_index)];
  sym_v2 = list->mem_fund.sym_index[vertex_to_index_fund(triang->vertices[1], list->mem_fund.vert_to_index)];
  sym_v3 = list->mem_fund.sym_index[vertex_to_index_fund(triang->vertices[2], list->mem_fund.vert_to_index)];
  for (int s=0; s < 48; s++) {
    tri_index index;
    if (triangle_in_fund(sym_v1[s], sym_v2[s], sym_v3[s], list->mem_fund.fund_len)) {//Point in the fund domain
      indices_unique_fund(sym_v1[s], sym_v2[s], sym_v3[s],index); 
      SMI(list->t_arr, index);
    } 
  } 
}

/*
 * Returns whether the triangle with vertices v1,v2 and v3 is inside the mem_list of type FUND.
 */
int mem_list_get_fund(tri_mem_list * list, arr3 v1, arr3 v2, arr3 v3){
  tri_index cur_index;
  
  //Indices of the vertices stored in cur_index
  cur_index[0] = vertex_to_index_fund(v1,list->mem_fund.vert_to_index);
  cur_index[1] = vertex_to_index_fund(v2,list->mem_fund.vert_to_index);
  cur_index[2] = vertex_to_index_fund(v3,list->mem_fund.vert_to_index);
  
  //Symmetry number needed to transform first vertex into fund domain
  int sym_num = list->mem_fund.vert_fund_sym[cur_index[0]]; 

  tri_index sym_index;
  
  //Apply symmetry_number and return the unique_index of this transformed triangle.
  indices_unique_fund(list->mem_fund.sym_index[cur_index[0]][sym_num], 
                      list->mem_fund.sym_index[cur_index[1]][sym_num],
                      list->mem_fund.sym_index[cur_index[2]][sym_num],
                      sym_index);
                        
  return GMI(list->t_arr,sym_index);
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
 * Converts an list of indices to a mem_list of cube type.
 * Sets all the symmetries of the indices as well.
 */
tri_mem_list mem_list_from_index_list(tri_index_list * list) {
  tri_mem_list result = mem_list_init_cube(list->dim, MEM_LIST_FALSE);
  //Now lets set all the nodes, exciting!
  for (size_t i = 0; i < list->len; i++) {
    triangle cur_triang = triangle_from_index_cube(list->index_list[i], list->dim);
    mem_list_set_cube(&result, &cur_triang);
  }
  return result;  
}

int index_list_from_file(tri_index_list * result, char * filename){
  FILE * stream;
  stream = fopen(filename, "rb");
  if (stream == NULL)
    return 0;
  if (fread(&result->dim, sizeof(int), 1, stream) < 1)
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
  fwrite(&list->dim, sizeof(int), 1, stream);
  fwrite(&list->len, sizeof(size_t), 1, stream);
  fwrite(list->index_list, sizeof(tri_index), list->len, stream);
  fclose(stream);
  return 1;  
}

