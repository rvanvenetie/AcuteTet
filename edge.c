#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <time.h>

#include "omp.h"
#include "triangle.h"
#include "edge.h"
#include "edge_list.h"


void print_edge(pedge edge) {
  printf("[[%d,%d],\n",edge->vertices[0][0],edge->vertices[0][1]);
  printf(" [%d,%d]]\n",edge->vertices[1][0],edge->vertices[1][1]);
}

//Returns whether this edge lies on the boundary of a [0..p] \times [0..p] 
#define mat2_col_equal(mat,col,val) (mat[0][col] == val && mat[1][col] == val)
int edge_boundary(pedge edge, int p) {
  return (mat2_col_equal(edge->vertices,0,0) || mat2_col_equal(edge->vertices,0,p) ||
          mat2_col_equal(edge->vertices,1,0) || mat2_col_equal(edge->vertices,1,p));
}

int edge_equal(pedge e1, pedge e2) {
	//AB = AB or AB = BA
	return ((equalArr2(e1->vertices[0], e2->vertices[0]) && equalArr2(e1->vertices[1], e2->vertices[1])) ||
			    (equalArr2(e1->vertices[1], e2->vertices[0]) && equalArr2(e1->vertices[0], e2->vertices[1])));
}




int edge_conform(pedge edge, edge_conform_parameters * parameters) {
	edge_matrix * mat = parameters->mat; //For ease of notation

	vert_index v1, v2; 
	v1 = vertex_to_index_square(edge->vertices[0], mat->p);
	v2 = vertex_to_index_square(edge->vertices[1], mat->p);
	if (!mat->val[v1][v2]) //This edge is not in the matrix
		return 0; 
	arr2 edge_normal; //Normal to this edge
  normalEdge(edge_normal, edge->vertices[0], edge->vertices[1]);
	int edge_d = dotArr2(edge_normal, edge->vertices[0]); //Find constant for this line
	
	int acute_above = 0, acute_below = 0;
	parameters->boundary_edge = edge_boundary(edge, mat->p);
	parameters->acute_above = 0;
	parameters->acute_below = 0;
	parameters->acute_ind_len = 0;
	
	//Loop over all points in the square
	for (vert_index i = 0; i < (mat->p + 1) * (mat->p + 1); i++) {
		arr2 point;
		vertex_from_index_square(point, i, mat->p);

		int dotprod = dotArr2(point, edge_normal); //Calculate <normal,point>, to find the side of this point
		if (((dotprod > edge_d && !acute_above) || //Point above
				(dotprod < edge_d && !acute_below))&&   //Point below
				triangle_acute_2d(edge->vertices[0], edge->vertices[1], point)  && //Acute
				mat->val[v1][i] && //Edge in mat
				mat->val[v2][i]) //Edge in mat
		{
			if (dotprod > edge_d)
				parameters->acute_above = 1;
			else
				parameters->acute_below = 1;

			if (parameters->store_acute_ind) //We want to store them all
				parameters->acute_ind[parameters->acute_ind_len++] = i;
			else { //Exit early if possible
				if ((parameters->acute_above && parameters->acute_below) || parameters->boundary_edge)
					return 1;

				acute_above = parameters->acute_above;
				acute_below = parameters->acute_below;
			}
		}
	}

	if (parameters->store_acute_ind) {

		if (parameters->acute_above && parameters->acute_below)
			return 1;

		if (parameters->boundary_edge && (parameters->acute_above || parameters->acute_below))
			return 1;
	}
	return 0;  
}

int edges_conform_sym_square(edge_matrix * mat) {
	int iterations = 0;
	int changed = 1;
	int row_size = edge_matrix_row_size(mat);
	double time_start =0 , time_end = 0;
	double time_entire = omp_get_wtime();
	square_points fund = gen_fund_square(mat->p);
	edge_conform_parameters parameters = {.mat = mat};
	while (changed) {
		iterations++;
		changed = 0;
		fprintf(stderr,"Starting conform loop with %zu edges.\n"
				, edge_matrix_count(mat));   
		time_start = omp_get_wtime();
		//Loop over fundamental domain
		for (int f = 0; f < fund.len; f++) {
			int i = vertex_to_index_square(fund.points[f], mat->p);
			for (int j = 0; j < row_size; j++)
				if (mat->val[i][j])
				{
				  struct edge	e;
					vertex_from_index_square(e.vertices[0], i, mat->p);
					vertex_from_index_square(e.vertices[1], j, mat->p);
					if (!edge_conform(&e, &parameters))
					{
						changed = 1;
						edge_matrix_clear_sym(mat, &e);
					} 
				}
		}
		time_end   = omp_get_wtime();
		fprintf(stderr,"Loop took %f seconds.\n\n", time_end-time_start);
	}
	fprintf(stderr,"Conforming set took %f seconds.\n", omp_get_wtime() - time_entire);
}

int edges_conform_square(edge_matrix * mat) {
	int iterations = 0;
	int changed = 1;
	int row_size = edge_matrix_row_size(mat);
	double time_start =0 , time_end = 0;
	double time_entire = omp_get_wtime();
	edge_conform_parameters parameters = {.mat = mat};
	while (changed) {
		iterations++;
		changed = 0;
		fprintf(stderr,"Starting conform loop with %zu edges.\n"
				, edge_matrix_count(mat));   
		time_start = omp_get_wtime();
		//Loop over upper triangular part (all possible edges)
		for (int i = 0; i < row_size; i++) 
			for (int j = i; j < row_size; j++)
				if (mat->val[i][j])
				{
				  struct edge	e;
					vertex_from_index_square(e.vertices[0], i, mat->p);
					vertex_from_index_square(e.vertices[1], j, mat->p);
					if (!edge_conform(&e, &parameters))
					{
						changed = 1;
						edge_matrix_clear(mat, &e);
					} 
				}
		time_end   = omp_get_wtime();
		fprintf(stderr,"Loop took %f seconds.\n\n", time_end-time_start);
	}
	fprintf(stderr,"Conforming set took %f seconds.\n", omp_get_wtime() - time_entire);
	return iterations;
}

#define TOUCH 2
#define DISJOINT 1
#define INTERSECT 0

int vert_vert_separation_axis(arr2 * vert1, int len1, arr2 * vert2, int len2, arr2 normal) {
  int min_1 = INT_MAX, max_1 = INT_MIN;
  int min_2 = INT_MAX, max_2 = INT_MIN;
	//Calculate projections of edge on normal
	for (int i = 0; i < len1; i++) {
		int dot1 = dotArr2(normal, vert1[i]);
		if (dot1 < min_1) min_1 = dot1;
		if (dot1 > max_1) max_1 = dot1;
	}
	//Calculate projections of triangle on normal
	for (int i = 0; i < len2; i++) {
		int dot2 = dotArr2(normal, vert2[i]);
		if (dot2 < min_2) min_2 = dot2;
		if (dot2 > max_2) max_2 = dot2;
	}

  /*
   * Disjoint if either range [min_1, max_1] is before [min_2, max_2] or
   * [min_2, max_2] is before [min_1, max_1]
   */
  if (max_1 < min_2 || max_2 < min_1)
    return DISJOINT;
  else if (max_1 > min_2 && max_2 > min_1) //!(max_1 <= min_2 || max_2 <= min_1)
    return INTERSECT;
  else { 
		//Both objects have vertices in the same line (touch).

		//Calculate the dotproduct of the line where we have touches
    int dot_line = (max_1 == min_2)? max_1 : max_2;

		//Count the amount of vertices on this line for edge/triangle
		int cnt_e = 0, cnt_t = 0;
		for (int i = 0; i < len1; i++)
			if (dotArr2(normal, vert1[i]) == dot_line)
				cnt_e++;
		for (int i = 0; i < len2; i++)
			if (dotArr2(normal, vert2[i]) == dot_line)
				cnt_t++;

		int cnt_min = (cnt_e < cnt_t)? cnt_e : cnt_t;
		
		//We should share  min_cnt vertices if we have a proper touch
		int cnt_share = 0;
		for (int i = 0; i < len1; i++) 
			for (int j = 0; j < len2; j++)
				if (equalArr2(vert1[i], vert2[j]))
					cnt_share++;

		if (cnt_share == cnt_min)
			return DISJOINT;
		else
      return INTERSECT;//TOUCH
  }
}
int edge_tri_separation_axis(pedge edge, ptriangle_2d triang, arr2 normal) {
	return vert_vert_separation_axis(edge->vertices, 2, triang->vertices,3,normal);
}

int tri_tri_separation_axis(ptriangle_2d t1, ptriangle_2d t2, arr2 normal) {
	return vert_vert_separation_axis(t1->vertices, 3, t2->vertices,3,normal);
}

int edge_tri_disjoint(pedge edge, ptriangle_2d triang) {
	arr2 edge_normal[3];
	//Check triangle for seperating axis
	normalEdge(edge_normal[0], triang->vertices[0], triang->vertices[1]);
	normalEdge(edge_normal[1], triang->vertices[0], triang->vertices[2]);
	normalEdge(edge_normal[2], triang->vertices[1], triang->vertices[2]);
	for (int i = 0; i < 3; i++)
		if (edge_tri_separation_axis(edge,triang,edge_normal[i]) == DISJOINT)
			return DISJOINT;

	//Check line segment (edge) for seperating axis
	normalEdge(edge_normal[0],edge->vertices[0], edge->vertices[1]);
	if (edge_tri_separation_axis(edge,triang,edge_normal[0]) == DISJOINT)
		return DISJOINT;

	return INTERSECT;
}

int tri_tri_disjoint(ptriangle_2d t1, ptriangle_2d t2) {
	arr2 edge_normal[3];
	//Check edges t1
	normalEdge(edge_normal[0], t1->vertices[0], t1->vertices[1]);
	normalEdge(edge_normal[1], t1->vertices[0], t1->vertices[2]);
	normalEdge(edge_normal[2], t1->vertices[1], t1->vertices[2]);
	for (int i = 0; i < 3; i++)
		if (tri_tri_separation_axis(t1,t2,edge_normal[i]) == DISJOINT)
			return DISJOINT;

	//Check edges t2
	normalEdge(edge_normal[0], t2->vertices[0], t2->vertices[1]);
	normalEdge(edge_normal[1], t2->vertices[0], t2->vertices[2]);
	normalEdge(edge_normal[2], t2->vertices[1], t2->vertices[2]);
	for (int i = 0; i < 3; i++)
		if (tri_tri_separation_axis(t1,t2,edge_normal[i]) == DISJOINT)
			return DISJOINT;

	return INTERSECT;
}


size_t filter_intersection_edge_matrix_tri(edge_matrix * mat, ptriangle_2d triang) {
	size_t removed = 0;
	int row_size = edge_matrix_row_size(mat);
	for (int i = 0; i < row_size; i++) 
		for (int j = i; j < row_size; j++)
			if (mat->val[i][j])
			{
				struct edge	e;
				vertex_from_index_square(e.vertices[0], i, mat->p);
				vertex_from_index_square(e.vertices[1], j, mat->p);
				if (!edge_tri_disjoint(&e, triang)) {
					removed++;
					edge_matrix_clear(mat, &e);
				}
			}
	return removed;
}

void print_triangulation(ptriangulation result) {
	printf("Triangulation debug:\n"
			   "\tp: %d\n"
				 "\tbound_len : %d\n"
				 "\tri_len : %d\n", result->p, result->bound_len, result->tri_len);
	printf("\nBoundaries:\n");
	for (int i = 0; i < result->bound_len; i++)
		print_edge(result->bound_edge + i);
	printf("\nTriangles:\n");
	for (int i = 0; i < result->tri_len; i++)
		print_triangle_2d(result->tri + i);
}

void free_triangulation(ptriangulation result) {
	free(result->bound_edge);
	free(result->tri);
}

//Returns whether the given triangle is disjoint with the given triangulation
int tri_triangulation_disjoint( ptriangle_2d tri, ptriangulation triang) {
	for (int i = 0; i < triang->tri_len; i++)
		if (!tri_tri_disjoint(triang->tri + i, tri))
			return INTERSECT;

	return DISJOINT;
}

//Removes all the triangles in the triangle list that intersect with triangulation, returns the new length
int filter_tri_list_intersect_triangulation(ptriangulation triang, ptriangle_2d tri_list, int tri_len) {
	int cntr = 0;
	for (int i = 0; i < tri_len; i++)
		if (tri_triangulation_disjoint(tri_list + i, triang))
			tri_list[cntr++] = tri_list[i];

	return cntr;
}
void add_boundary_triangulation(ptriangulation result, arr2 v1, arr2 v2) {
	struct edge e = (struct edge) {{{v1[0], v1[1]}, {v2[0], v2[1]}}};
	//Check if this boundary is already contained in our boundary array
	//If so, we have to remove this, as we do not need further expansion here
	for (int i = 0; i < result->bound_len; i++)
		if (edge_equal(&e, result->bound_edge + i)) {
			memmove(result->bound_edge + i, result->bound_edge + i + 1, (--result->bound_len - i) * sizeof(struct edge));
			return;
		}

	if (edge_boundary(&e, result->p)) //Edge lies on boundary square, we do not add this one
		return;
	result->bound_edge = realloc(result->bound_edge, (++result->bound_len) * sizeof(struct edge));
	result->bound_edge[ result->bound_len - 1] = e;
}

void add_tri_triangulation(ptriangulation result, ptriangle_2d triang) {
	result->tri = realloc(result->tri, (++result->tri_len) * sizeof(triangle_2d));
	result->tri[result->tri_len - 1] = *triang;

	//Add edges to the boundary
	add_boundary_triangulation(result, triang->vertices[0], triang->vertices[1]);
	add_boundary_triangulation(result, triang->vertices[1], triang->vertices[2]);
	add_boundary_triangulation(result, triang->vertices[0], triang->vertices[2]);
} 

triangulation triangulate_square(edge_matrix * mat) {
	triangulation result = {0};
	result.p = mat->p;

	//Find a starting edge on the boundary (0,0)..(x, 0)
	result.bound_edge = calloc(sizeof(struct edge), 1);
	do {
		result.bound_edge->vertices[1][0] = rand() % (mat->p+1);
	} while (! edge_matrix_contains(mat, result.bound_edge));
	result.bound_len = 1;

	//Initalization
	edge_conform_parameters parameters = {
		.mat = mat, 
		.store_acute_ind = 1, 
		.acute_ind = malloc(sizeof(vert_index) * (mat->p + 1) * (mat->p + 1))
	};
	
  fprintf(stderr,"Starting triangulation with edge:\n");
	//print_edge(result.bound_edge);

	//While edges on the boundary
	while (result.bound_len > 0) {
    /*
     * We are going to add a tetrahedron on the boundary triangle.
     * To do so, we select a random triangle on the boundary. Then we generate all the
     * acute tetrahedra (above and below) with facets in our possible list.
     * From this list we remove all the tetrahedrons that intersect with our current triangulation.
     * Then we add a random tetrahedron to our triangulation, update the conform list and repeat.
     */
    int rand_bound = rand() % result.bound_len;
		fprintf(stderr,"\n\nTotal amount of edges left:%zu\n", edge_matrix_count(mat));
		fprintf(stderr,"Expanding triangulation at boundary edge :\n");
		//print_edge(result.bound_edge + rand_bound);

		if (!edge_conform(&result.bound_edge[rand_bound], &parameters)) 
    {
			//Happens by conforming the edge_matrix
      fprintf(stderr,"We have an edge on the boundary that is not conform anymore.\n");
      break;
    }
		
		int 			tri_list_len;
		ptriangle_2d tri_list;

		tri_list_len = parameters.acute_ind_len;
    fprintf(stderr,"Total amount of conform triangles found for this boundary: %hu\n", tri_list_len);

		//Form explicit list of triangles
		tri_list = malloc(sizeof(triangle_2d) * tri_list_len);
    for (unsigned short i = 0; i < tri_list_len; i++) 
    {
			copyArr2(tri_list[i].vertices[0], result.bound_edge[rand_bound].vertices[0]);
			copyArr2(tri_list[i].vertices[1], result.bound_edge[rand_bound].vertices[1]);
			vertex_from_index_square(tri_list[i].vertices[2], parameters.acute_ind[i], mat->p);
    }

		/*
		 * Remove all triangles that intersect with the current triangulation
		 */
    tri_list_len = filter_tri_list_intersect_triangulation(&result, tri_list, tri_list_len);
		if (!tri_list_len) 
		{
			//Happens, triangulation gets stuck somewhere..
			fprintf(stderr,"Dafuq!? All triangles intersect with our triangulation..\n");
			break;
		}
		//Select random triangle from acute list
		int rand_tri = rand() % tri_list_len;
		fprintf(stderr,"Adding following triangle to the triangulation\n");
		//print_triangle_2d(tri_list + rand_tri);
		add_tri_triangulation( &result,tri_list + rand_tri);
    if (!result.bound_len) //If we have no boundaries left, we must be done!!
    {
      fprintf(stderr,"No more boundaries left.. WE FINNISHED!??\n");
      break;
    }

		//Remove all edges interescting with this triangle
		size_t removed = filter_intersection_edge_matrix_tri(mat, tri_list + rand_tri);
		fprintf(stderr,"Removed %zu edges intersecting with the new triangle\n", removed);

		//Conform edge matrix
		edges_conform_square(mat);
  }
	//print_triangulation(&result);
  free(parameters.acute_ind);
  fprintf(stderr,"Triangulation has length of %zu\n", result.tri_len);
  return result;
}

#define get_cpu_time(start) (((double) (clock() - start)) / CLOCKS_PER_SEC)
int main(int argc, char *argv[]) {
	edge_matrix mat, mat_copy;
	clock_t start;
	double time_init, time_conform, time_conform_sym, time_cosy, time_triang, suc_rate;
	size_t edge_total, edge_con_cnt, cosy_tri, cosy_area;
	int p = 16;
	printf("%-5s" "%-15s"       "%-15s"					"%-15s"       "%-15s"        "%-15s"			"%-15s"			"%-15s"		"%-15s"    "%-15s" "\n",
			   "p",   "edge_total","edge_conform","time_conf","time_conf_sym","cosy_tri","cosy_area_x2","time_cosy","suc_rate","time_triang");
	//printf("p\tedge_total\tedge_con_cnt\tcosy_tri\tcosy_area\n");
	for (int p = 1; p < 40; p++)
	{
		start = clock();
		mat = edge_matrix_init(p, 1);
		time_init = get_cpu_time(start);

		mat_copy = edge_matrix_copy(&mat);
		edge_total = edge_matrix_count(&mat);

		start = clock();
		edges_conform_square(&mat);
		time_conform = get_cpu_time(start);

		start = clock();
		edges_conform_sym_square(&mat_copy);
		time_conform_sym = get_cpu_time(start);

		edge_con_cnt = edge_matrix_count(&mat);
		if (edge_con_cnt != edge_matrix_count(&mat_copy))
		{
			printf("ERROR ERROR, SYMMETRY FINDS OTHER COUNT\n");
			printf("%zu-%zu\n", edge_con_cnt, edge_matrix_count(&mat_copy));
		}

		start = clock();
		edge_matrix_cosy_count(&mat, &cosy_tri, &cosy_area);
		time_cosy = get_cpu_time(start);

		edge_matrix_free(&mat_copy);
		triangulation triang;
		int suc_cntr = 0, total_cntr = 0;
		start = clock();
    while (suc_cntr != 20 && edge_con_cnt > 0) {
			total_cntr++;
			mat_copy = edge_matrix_copy(&mat);
			triang = triangulate_square(&mat_copy);
			edge_matrix_free(&mat_copy);
			if (triang.bound_len == 0) //Succesful triangulation
				suc_cntr++;
			free_triangulation(&triang);
		}
		time_triang = get_cpu_time(start);
		 suc_rate = (double) suc_cntr / total_cntr;
		//print_triangulation(&triang);
		edge_matrix_free(&mat);

		fprintf(stderr,"\n\nCalculating edge matrix for p = %d\n", p);
		fprintf(stderr,"Edge count total  = %zu\n", edge_total);
		fprintf(stderr,"Edge count mat    = %zu\n", edge_con_cnt);
		fprintf(stderr,"Total cosy tri    = %zu\n", cosy_tri);
		fprintf(stderr,"Total cosy area   = %zu.%d\n", cosy_area/2,5 * (cosy_area % 2));
		fprintf(stderr,"Succes rate triang= %f\n", suc_rate);
	   printf("%-5d" "%-15zu" "%-15zu"      "%-15f"      "%-15f"          "%-15zu"  "%-15zu"   "%-15f"   "%-15f"   "%-15f" "\n",
				    p,   edge_total, edge_con_cnt,time_conform,time_conform_sym,cosy_tri, cosy_area, time_cosy,suc_rate, time_triang);
		//printf("%d\t%zu\t%zu\t%zu\t%zu.%d\n", p, edge_total, edge_con_cnt, cosy_tri, cosy_area/2, 5 * (cosy_area % 2));
	}
	return 0;
}
