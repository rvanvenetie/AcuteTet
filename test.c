#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "vector.h"
#include "triangle.h"
#include "tetraeder.h"
#include "mem_list.h"
#include "triangulate.h"
#include "combinations.h"

void test_sym(void){
  for (int i = 1; i < 100; i++) {
    int dim = i;
    size_t dim_size = (i+1) * (i+1) * (i+1);
    unsigned char * points = calloc(dim_size, sizeof(unsigned char)); 
    cube_points fund = gen_fund_points(i);
    for (size_t j = 0; j < fund.len; j ++) //For each point in the fundamental domain
      for (int k = 0; k < 48; k++) {//Apply symmetry
        arr3 sym_fund;
        apply_symmetry(k,i, fund.points[j], sym_fund);
        points[vertex_to_index_cube(sym_fund, dim)] = 1;
      }
    for (size_t j = 0; j < dim_size; j++) {
      if (points[j] == 0) {
        printf("Symmetry test failed. Dimension %d.\n",i);
        exit(0);
      }
    }
    free(points);
  }
}

triangle rand_triangle(int dim) {
  triangle result;
  dim = dim + 1;
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

void rand_vertex_tet(int dim, arr3 vertex) {
  vertex[0] = rand() % dim;
  vertex[1] = rand() % (dim - vertex[0]);
  vertex[2] = rand() % (dim - vertex[0] - vertex[1]);
}
triangle rand_triangle_tet(int dim) {
  triangle result;
  rand_vertex_tet(dim, result.vertices[0]);
  rand_vertex_tet(dim, result.vertices[1]);
  rand_vertex_tet(dim, result.vertices[2]);
  if (equalArr3(result.vertices[0], result.vertices[1]) ||
      equalArr3(result.vertices[1], result.vertices[2]) ||
      equalArr3(result.vertices[0], result.vertices[2]))
    return rand_triangle_tet(dim);
  return result;
}


tetra rand_tetra(int dim) {
  tetra result;
  result.vertices[0][0] = rand() % dim;
  result.vertices[0][1] = rand() % dim;
  result.vertices[0][2] = rand() % dim;
  result.vertices[1][0] = rand() % dim;
  result.vertices[1][1] = rand() % dim;
  result.vertices[1][2] = rand() % dim;
  result.vertices[2][0] = rand() % dim;
  result.vertices[2][1] = rand() % dim;
  result.vertices[2][2] = rand() % dim;
  result.vertices[3][0] = rand() % dim;
  result.vertices[3][1] = rand() % dim;
  result.vertices[3][2] = rand() % dim;
  return result;
}

int tri_index_same(tri_index a, tri_index b){
  return (a[0] == b[0] && a[1] == b[1] && a[2] == b[2]);
}
void test_triangle_indices(void){
  int dim = 50;
  for (int i = 0; i < 500; i++){
    triangle origin = rand_triangle(dim);
    tri_index ori,v1,v2;
    triangle_to_index_cube(origin, dim, ori);
    vertices_to_index_cube(origin.vertices[1], origin.vertices[2], origin.vertices[0], dim, v1);
    vertices_to_index_cube(origin.vertices[1], origin.vertices[0], origin.vertices[2], dim, v2);
    if (!tri_index_same(ori,v1) || !tri_index_same(v1,v2)) {
      printf("UNIQUE INDEX NOT THE SAME!!\n");
      print_triangle(&origin);
    }
  }
}
void swap_vertex(arr3 a, arr3 b) {
  arr3 t;
  t[0] = b[0];
  t[1] = b[1];
  t[2] = b[2];
  
  b[0] = a[0];
  b[1] = a[1];
  b[2] = a[2];
  
  a[0] = t[0];
  a[1] = t[1];
  a[2] = t[2];
}
void randomize_triangle(ptriangle triang) {
  swap_vertex(triang->vertices[rand() % 3], triang->vertices[rand() % 3]);
  swap_vertex(triang->vertices[rand() % 3], triang->vertices[rand() % 3]);
  swap_vertex(triang->vertices[rand() % 3], triang->vertices[rand() % 3]);
  swap_vertex(triang->vertices[rand() % 3], triang->vertices[rand() % 3]);
  swap_vertex(triang->vertices[rand() % 3], triang->vertices[rand() % 3]);
}
void test_mem_list_fund(void){
  int dim = 15;
  tri_mem_list mem_list = mem_list_init_fund(dim,MEM_LIST_FALSE);
  triangle * triang_list = malloc(sizeof(triangle)  * 750);
  for (int i = 0; i < 750; i++) {
    triang_list[i] = rand_triangle(dim);
    mem_list_set_fund(&mem_list,triang_list + i);
  } //Inserted all the triangles.
  for (int i = 0; i < 750; i++) {
    randomize_triangle(triang_list + i);
    triangle_symmetry(triang_list + i, rand() % 48,dim, triang_list + i); 
    if (!mem_list_get_fund(&mem_list, triang_list[i].vertices[2], triang_list[i].vertices[1], triang_list[i].vertices[0])) {
      printf("Triangle not fount!!");
    }
  } //Check for additional 
  free(triang_list);
}

void test_mem_list_tet(void) {
  int dim = 15;
  int len = 2500;
  tri_mem_list mem_list = mem_list_init_tet(dim,MEM_LIST_FALSE);
  triangle * triang_list = malloc(sizeof(triangle)  * len);
  for (int i = 0; i < len; i++) {
    triang_list[i] = rand_triangle_tet(dim);
    mem_list_set_tet(&mem_list,triang_list + i);
  } //Inserted all the triangles
  printf("Amount of unique triangles: %zu\n", mem_list_count(&mem_list));
  for (int i = 0; i < len; i++) {
    randomize_triangle(triang_list + i);
    if (!mem_list_get_tet(&mem_list, triang_list[i].vertices[2], triang_list[i].vertices[1], triang_list[i].vertices[0])) {
      printf("Triangle not fount!!");
    }
  } 
  for (int i = 0; i < len; i++)  {
    mem_list_clear_tet(&mem_list, triang_list + i);
  }
  printf("Amount of triangles left : %zu\n", mem_list_count(&mem_list));
  mem_list_free(&mem_list);
  free(triang_list);
}

int tetra_acute_normal(ptetra tet) {
  arr3 P[5]; //Edges
  arr3 normals[4]; 
  /*
   * Calculate three edges of the tetrahedron
   */
  subArr3(tet->vertices[1], tet->vertices[0], P[0]);
  subArr3(tet->vertices[2], tet->vertices[0], P[1]);
  subArr3(tet->vertices[3], tet->vertices[0], P[2]);
  
  /*
   * All Facets must have acute triangles, cheap test to rule this tetra out right now
   */
  if (dotArr3(P[0],P[1]) < 0 || dotArr3(P[0], P[2]) < 0 || dotArr3(P[1],P[2]) < 0 )
    return 0;
  
  crossArr3(P[2],P[0], normals[2]); //Normal on facet 0,1,3
  crossArr3(P[0],P[1], normals[3]); //Normal on facet 0,1,2
  
  if (dotArr3(normals[2], normals[3]) >= 0)
    return 0;
    
  crossArr3(P[1],P[2], normals[1]); //Normal on facet 0,2,3
  if (dotArr3(normals[1], normals[2]) >= 0 ||
      dotArr3(normals[1], normals[3]) >= 0)
    return 0;
  
  subArr3(tet->vertices[2], tet->vertices[1], P[3]);
  subArr3(tet->vertices[3], tet->vertices[1], P[4]); 
  crossArr3(P[4],P[3], normals[0]); //Normal on facet 1,2,3 
  
  return (dotArr3(normals[0], normals[1]) < 0 &&
          dotArr3(normals[0], normals[2]) < 0 &&
          dotArr3(normals[0], normals[3]) < 0);
}

void test_tetra_normals(void) {
  arr3 normals[4];
  for (int i = 0; i < 1000000; i++) {
    tetra tet = rand_tetra(20); 
    
    tetra_normals(&tet, normals);
    int pos = 0;
    int neg = 0;
    for (int j = 0; j < 4; j++) { //Test orientation of normal (check against apex)
      arr3 test_vec;
      /*
       * Normals[j] is the facet without vertex j.
       * Thus the apex is vertex [j], and point in the facet is given by any other point than j
       */
      subArr3(tet.vertices[j], tet.vertices[(j+1)%4], test_vec); //Apex - face[j]
      if (dotArr3(test_vec, normals[j]) > 0)
        pos = 1;
      else
        neg = 1;
    }
    if (pos && neg) {
      printf("Inconsistent normals on tetra!\n");
      print_tetra(&tet);
    }
    if (tetra_acute(&tet) ^ tetra_acute_normal(&tet)) //XOR, means we got different results
    {
      printf("Tetra acute by one method, but not by the other.\n");
      printf("Tetra_acute method: %d\n", tetra_acute(&tet));
      printf("Test acute method : %d\n", tetra_acute_normal(&tet));
      print_tetra(&tet);
    }
    
  }
}

void test_tetra_disjunct(void) {
  tetra t1 = (tetra) {{{5,5,5},{6,5,5},{5,6,5},{5,5,6}}};
  tetra t2 = (tetra) {{{5,5,5},{6,5,5},{5,6,5},{5,5,4}}};
  printf("Test disjunct: %d", tetra_tetra_disjoint(&t1,&t2));
}

void test_boundary(int dim) {
  for (int i = 0; i < 15000; i ++) {
    triangle tmp_tri = rand_triangle(dim);
    if (triangle_boundary_cube(&tmp_tri,dim))
      print_triangle(&tmp_tri);
  }
}

void test_prism(int dim) {
  while (1){
    tetra tmp_tetra = rand_tetra(dim);
    mat3 P;
    triangle_sides(tmp_tetra.vertices[0], tmp_tetra.vertices[1], tmp_tetra.vertices[2],P);
    arr3 tri_normal;
    crossArr3(P[1], P[0], tri_normal); //Calculate normal on the triangle plane
    //Calculate the vector perpendicular to each side and the normal. Normals on each side of the prism
    mat3 side_normals;
    arr3 side_d; //The constant expression for the plane equations of the sides
    //Convention, need to explain why it works?? Third point must be on the other side!!
  
    crossArr3(tri_normal, P[0], side_normals[0]);
    crossArr3(P[1], tri_normal, side_normals[1]);
    crossArr3(tri_normal, P[2], side_normals[2]);
    side_d[0] = dotArr3(side_normals[0], tmp_tetra.vertices[0]); 
    side_d[1] = dotArr3(side_normals[1], tmp_tetra.vertices[0]);
    side_d[2] = dotArr3(side_normals[2], tmp_tetra.vertices[2]);
    int a = (dotArr3(tmp_tetra.vertices[3],side_normals[0]) < side_d[0] &&
          dotArr3(tmp_tetra.vertices[3],side_normals[1]) < side_d[1] &&
          dotArr3(tmp_tetra.vertices[3],side_normals[2]) < side_d[2]);
    int b = tetra_acute(&tmp_tetra);
    if (b && !a)
    printf("%d,%d\n",a,b);
  }
}
int main(void){
  test_prism(20);
  //test_boundary(5);
  //test_sym();
  //test_tetra_normals();
  //test_mem_list_tet();
  //test_triangle_indices();
  //
  //test_tetra_disjunct();
}
