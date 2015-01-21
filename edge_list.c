#include <stdio.h>
#include <string.h>

#include "triangle.h"
#include "vector.h"
#include "edge.h"
#include "edge_list.h"

#ifndef INLINE_MACROS
vert_index vertex_to_index_square(arr2 vert, int p) {
	return (vert[0] * (p + 1) + vert[1]);
}

vert_index xy_to_index_square(int x, int y, int p) {
	return (x * (p+1) + y);
}

void  vertex_from_index_square(arr2 result, vert_index index, int p) {
	result[0] = index / (p + 1);
	result[1] = index % (p + 1);
}

void edge_to_index_square(edge_index result, pedge edge, int p) {
	result[0] = vertex_to_index_square(edge->vertices[0], p);
	result[1] = vertex_to_index_square(edge->vertices[1], p);
}

void edge_from_index_square(pedge result, edge_index indices, int p) {
	vertex_from_index_square(result->vertices[0], indices[0], p);
	vertex_from_index_square(result->vertices[1], indices[1], p);
}
#endif
/*
 * Returns the row size for this edge matrix. (Amount of possible edges)
 */
int edge_matrix_row_size(edge_matrix * mat) {
	return ((mat->p + 1) * (mat->p + 1));
}

#define EDGE_INIT_TRUE 1
#define EDGE_INIT_FALSE 0
edge_matrix edge_matrix_init(int p, unsigned char init_val) {
	edge_matrix result;
	result.p = p;
  int row_size = edge_matrix_row_size(&result);
	result.val = malloc(sizeof(unsigned char *)* row_size);
  for (int i = 0; i < row_size; i++) {
		result.val[i] = malloc(sizeof(unsigned char) * row_size);
		memset(result.val[i], init_val, row_size);
		result.val[i][i] = 0;
	}
	return result;
}

edge_matrix edge_matrix_copy(edge_matrix * mat) {
	edge_matrix result;
	result.p = mat->p;
  int row_size = edge_matrix_row_size(&result);
	result.val = malloc(sizeof(unsigned char *)* row_size);
  for (int i = 0; i < row_size; i++) {
		result.val[i] = malloc(sizeof(unsigned char) * row_size);
		memcpy(result.val[i], mat->val[i], row_size);
	}
	return result;

}
void edge_matrix_free(edge_matrix * mat) {
  int row_size = edge_matrix_row_size(mat);
	for (int i  =0; i < row_size; i++)
		free(mat->val[i]);
	free(mat->val);
}

void edge_matrix_set(edge_matrix * mat, pedge edge) {
	edge_index indices;
	edge_to_index_square(indices,edge, mat->p);
	mat->val[indices[0]][indices[1]] = mat->val[indices[1]][indices[0]] = 1;
}

void edge_matrix_clear(edge_matrix * mat, pedge edge) {
	edge_index indices;
	edge_to_index_square(indices,edge, mat->p);
	mat->val[indices[0]][indices[1]] = mat->val[indices[1]][indices[0]] = 0;
}

void edge_matrix_clear_sym(edge_matrix * mat, pedge edge) {
	int AX, AY, BX, BY,N; //For easy of notation
	AX = edge->vertices[0][0];
	AY = edge->vertices[0][1];
	BX = edge->vertices[1][0];
	BY = edge->vertices[1][1];
	N = mat->p;

	int symmetries[8][4] = {
			{AX  ,AY  ,BX  ,  BY}, //Original
			{AY  ,AX  ,BY  ,  BX}, //Mirror in the line y = x
			{AY  ,N-AX,BY  ,N-BX}, //Flip vertical y' = N - y, x' = x
			{AX  ,N-AY,BX  ,N-BY}, //Mirror in line y' = N -x, x' = N - y
			{N-AX,N-AY,N-BX,N-BY}, //Flip horizontal, y' = y, x' = M - x
			{N-AY,N-AX,N-BY,N-BX}, //Mirror in the line y= x
			{N-AY,AX  ,N-BY,  BX}, //Mirror vertical y' = N - y
			{N-AX,AY  ,N-BX,  BY}  //Mirror in the line y' = N - x
	};
	edge_index indices;
	for (int sym = 0; sym < 8; sym++) {
		indices[0] = xy_to_index_square(symmetries[sym][0], symmetries[sym][1], mat->p);
		indices[1] = xy_to_index_square(symmetries[sym][2], symmetries[sym][3], mat->p);
		mat->val[indices[0]][indices[1]] = mat->val[indices[1]][indices[0]] = 0;
	}
}

int  edge_matrix_contains(edge_matrix * mat, pedge edge) {
	edge_index indices;
	edge_to_index_square(indices,edge, mat->p);
	return (mat->val[indices[0]][indices[1]]);
}

unsigned char edge_matrix_get(edge_matrix * mat, arr2 v1, arr2 v2) {
	edge_index indices;
	indices[0] = vertex_to_index_square(v1,mat->p);
	indices[1] = vertex_to_index_square(v2,mat->p);
	return (mat->val[indices[0]][indices[1]]);
}

size_t edge_matrix_count(edge_matrix * mat ) {
	size_t result = 0;
	int row_size = edge_matrix_row_size(mat);
	//Loop over upper triangular part
	for (int i = 0; i < row_size; i++) 
		for (int j = i; j < row_size; j++)
			result += mat->val[i][j];
	return result;
}

void edge_matrix_cosy_count(edge_matrix * mat, size_t * tri_cnt, size_t * total_area_x2) {
	*tri_cnt = 0;
	*total_area_x2 = 0;
	//Amount of vertices in this matrix
	int row_size = edge_matrix_row_size(mat); 
	arr2 v1,v2,v3;
	//Loop over all triangles in the square
	for (int i = 0; i < row_size; i++) {
		vertex_from_index_square(v1, i, mat->p);
		for (int j = i; j < row_size; j++) {
			vertex_from_index_square(v2, j, mat->p);
			for (int k = j; k < row_size; k++) {
				vertex_from_index_square(v3, k, mat->p);
				//If triangle is acute and the sides are contained in the matrix
				if (triangle_acute_2d(v1,v2,v3) && 
						mat->val[i][j] &&
						mat->val[i][k] &&
						mat->val[j][k])
				{
					(*tri_cnt)++;
					(*total_area_x2) += triangle_area_x2(v1,v2,v3);
				}

			}
		}
	}
}

int edge_matrix_to_file(edge_matrix * mat, char * filename) {
  FILE * stream;
  stream = fopen(filename, "wb");
  if (stream == NULL)
    return 0;
  //Write p to the file
  if (fwrite(&mat->p, sizeof(int), 1, stream) < 1)
    return 0;
  
  //Save entire matrix to the file (possibly only upper triangular part?)
  int row_size = edge_matrix_row_size(mat);
	for (int i = 0; i < row_size; i++)  
    if (fwrite(mat->val[i], sizeof(unsigned char), row_size , stream) < (row_size ))
      return 0;

  fclose(stream);
  return 1;
}

int edge_matrix_from_file(edge_matrix * result, char  * filename) { 
  FILE * stream;
  stream = fopen(filename, "rb");
  if (stream == NULL)
    return 0;
  

  //Read p from the file
  if (fread(&result->p, sizeof(int), 1, stream) < 1)
    return 0;

  *result = edge_matrix_init(result->p, 0);
  //Read entire matrix from file
  int row_size = edge_matrix_row_size(result);
	for (int i = 0; i < row_size; i++)  
    if (fread(result->val[i], sizeof(unsigned char), row_size , stream) < (row_size ))
      return 0;
  fclose(stream);
  return 1;
}
