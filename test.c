#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "vector.h"
#include "triangle.h"
#include "tetraeder.h"
#include "mem_list.h"
#include "tri_list.h"
#include "triangulation.h"


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
  printf("Passed symmetry test\n");
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
  printf("Unique index cube-test passed\n");
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
  tri_mem_list mem_list = mem_list_fund_init(dim,MEM_LIST_FALSE);
  triangle * triang_list = malloc(sizeof(triangle)  * 750);
  for (int i = 0; i < 750; i++) {
    triang_list[i] = rand_triangle(dim);
    mem_list_fund_set(&mem_list,triang_list + i);
  } //Inserted all the triangles.
  
  //Check if all set
  for (int i = 0; i < 750; i++) {
    randomize_triangle(triang_list + i);
    triangle_symmetry(triang_list + i, rand() % 48,dim, triang_list + i); 
    if (!mem_list_fund_get(&mem_list, triang_list[i].vertices[2], triang_list[i].vertices[1], triang_list[i].vertices[0])) {
      printf("ERROR IN MEM_LIST_FUND");
    }
  } 
  //Delete all
  for (int i = 0; i < 750; i++) 
    mem_list_fund_clear(&mem_list, triang_list + i);
  //Check if all deleted
  if (mem_list_count(&mem_list) > 0)
    printf("ERROR 2 IN MEM_LIST_FUND\n");
  //Check for additional 
  free(triang_list);
  printf("MEM_LIST_FUND test passed\n");
}

void test_mem_list_tet(void) {
  int dim = 15;
  int len = 2500;
  tri_mem_list mem_list = mem_list_tet_init(dim,MEM_LIST_FALSE);
  triangle * triang_list = malloc(sizeof(triangle)  * len);
  for (int i = 0; i < len; i++) {
    triang_list[i] = rand_triangle_tet(dim);
    mem_list_tet_set(&mem_list,triang_list + i);
  } //Inserted all the triangles

  for (int i = 0; i < len; i++) {
    randomize_triangle(triang_list + i);
    if (!mem_list_tet_get(&mem_list, triang_list[i].vertices[2], triang_list[i].vertices[1], triang_list[i].vertices[0])) {
      printf("ERROR 1 IN MEM_LIST_TET");
    }
  } 
  for (int i = 0; i < len; i++)  {
    mem_list_tet_clear(&mem_list, triang_list + i);
  }
  if (mem_list_count(&mem_list) > 0)
    printf("ERROR 2 IN MEM_LIST_TET\n");
  mem_list_free(&mem_list);
  free(triang_list);
  printf("MEM_LIST_TET test passed\n");
}
long long C(int n, int r) {
  if(r > n / 2) r = n - r; // because C(n, r) == C(n, n - r)
  long long ans = 1;
  int i;

  for(i = 1; i <= r; i++) {
    ans *= n - r + i;
    ans /= i;
  }

  return ans;
}
void test_tri_list(void) {
  int dim = 3;
  int len = 2500;
  tri_list list = tri_list_init(dim, MEM_LIST_FALSE);
  if (tri_list_count(&list) != 0)
    printf("BLABLA ERROR TRI_LIST\n");
  triangle * triang_list = malloc(sizeof(triangle)  * len);
  for (int i = 0; i < len; i++) {
    triang_list[i] = rand_triangle_tet(dim);
    tri_list_insert(&list, triang_list + i,TRI_LIST_NO_RESIZE);
    triang_list[i] = rand_triangle_tet(dim);
    tri_list_insert(&list, triang_list + i,TRI_LIST_NO_RESIZE);
  } //Inserted all the triangles + some random other shit
  size_t mem = tri_list_memory(&list);
  tri_list_empty(&list); //Empty the list.. Should leave some space, so reinserting should not create any new memory.
  if (tri_list_count(&list))
    printf("ERROR TRI SHOULD BE EMPTY\n");
  for (int i = 0; i < len; i++)
    tri_list_insert(&list, triang_list + i, TRI_LIST_NO_RESIZE);

  if (tri_list_memory(&list) != mem)
    printf("MEMORY SHOULD NOT BE DIFFERENT");

  for (int i = 0; i < len; i++)
    if (!tri_list_contains(&list, triang_list +i))
      printf("ERROR, LIST SHOULD CONTAIN THIS TRIANGLE\n");
  size_t count = tri_list_count(&list);
  tri_list_to_file(&list, "test.tmp");
  tri_list_free(&list);
  tri_list_from_file(&list, "test.tmp");
  if (count != tri_list_count(&list))
    printf("ERROR 3 IN TRI_LIST");
  for (int i = 0; i < len; i++) {
    randomize_triangle(triang_list + i);
    if (!tri_list_contains(&list, triang_list + i)) {
      printf("ERROR 1 IN TRI_LIST");
    }
  } 
  for (int i = 0; i < len; i++)  {
    tri_list_remove(&list, triang_list + i, TRI_LIST_NO_RESIZE );
  }
  if (tri_list_count(&list) > 0)
    printf("ERROR 2 IN TRI_LIST\n");
  tri_list_free(&list);
  free(triang_list);

  dim = 5;
  tri_mem_list mem_list = mem_list_fund_init(dim, MEM_LIST_TRUE);
  list = mem_list_to_tri_list(&mem_list);
  if (C((dim + 1) * (dim + 1)  * (dim + 1), 3) != tri_list_count(&list))
    printf("ERROR 4 IN TRI_LIST\n");
  mem_list_free(&mem_list);
  tri_list_free(&list);

  mem_list = mem_list_fund_init(dim, MEM_LIST_FALSE);
  triangle start_triangle = rand_triangle_tet(dim);
  mem_list_fund_set(&mem_list, &start_triangle);
  printf("Mem_list_count: %zu\n", mem_list_count(&mem_list));
  list = mem_list_to_tri_list(&mem_list);
  printf("Tri_list count: %zu\n", tri_list_count(&list));

  for (int i = 0; i < 48; i++)
  {
    triangle sym_triangle;
    triangle_symmetry(&start_triangle, i, dim, &sym_triangle);
    if (!tri_list_contains(&list, &sym_triangle))
      printf("ERROR 5 IN TRI_LIST\n");
  }
  for (int i = 0; i < 48; i++)
  {
    triangle sym_triangle;
    triangle_symmetry(&start_triangle, i, dim, &sym_triangle);
    tri_list_remove(&list, &sym_triangle, TRI_LIST_NO_RESIZE);
  }

  if (tri_list_count(&list) != 0)
    printf("ERROR 6 IN TRI_LIST\n");
  printf("TRI_LIST test passed\n");
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
    
  }
  printf("Normals on tetrahedra-test passed\n");
}


void test_boundary(int dim) {
  for (int i = 0; i < 15000; i ++) {
    triangle tmp_tri = rand_triangle(dim);
    if (triangle_boundary_cube(&tmp_tri,dim)) {
      //print_triangle(&tmp_tri);
    }
  }
}

void test_prism(int dim) {
  for (int i = 0; i < 15000; i++){
    tetra tmp_tetra = rand_tetra(dim);
    arr3 P[3];
    triangle_sides(tmp_tetra.vertices[0], tmp_tetra.vertices[1], tmp_tetra.vertices[2],P);
    arr3 tri_normal;
    crossArr3(P[1], P[0], tri_normal); //Calculate normal on the triangle plane
    //Calculate the vector perpendicular to each side and the normal. Normals on each side of the prism
    arr3 side_normals[3];
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
    
    
    triangle tmp_triangle;
    memcpy(tmp_triangle.vertices, tmp_tetra.vertices, 3 * sizeof(arr3));
    int b = tetra_acute_optimized(&tmp_triangle, tmp_tetra.vertices[3]);
    if (b && !a)
    printf("ERROR IN PRISM TEST: %d,%d\n",a,b);
  }
  printf("Projection on facet-test passed\n");
}


void test_tetra_disjoint(void) {
  tetra tet;
  triangle tri;

  tet = (tetra) {{{0,0,0},{16,0,0},{12,9,0},{10,5,9}}};
  tri =(triangle) {{{0,0,0},{16,0,0},{12,9,0}}};
  //This triangle and tetrahedron share a facet. So they must be disjoint
  if (!tri_tet_disjoint(&tri, &tet))
    printf("ERROR1IN DISJOINT TEST");


  tri =(triangle) {{{0,0,0},{8,8,9},{12,9,0}}};
  //This triangle and tetrahedron share an edge and are disjoint
  if (!tri_tet_disjoint(&tri, &tet)) 
    printf("ERROR2IN DISJOINT TEST");

  tri =(triangle) {{{0,0,0},{12,5,9},{12,9,0}}};
  //This triangle and tetrahedron share an edge, but not disjoint, third point outside the tetrahedron
  if (tri_tet_disjoint(&tri, &tet)) 
    printf("ERROR3IN DISJOINT TEST");

  tri =(triangle) {{{0,0,0},{12,5,1},{12,9,0}}};
  //This triangle and tetrahedron share an edge, but not disjoint, third point inside the tetrahedron
  if (tri_tet_disjoint(&tri, &tet))
    printf("ERROR4IN DISJOINT TEST");

  tet = (tetra) {{{0,0,0},{2,0,0},{1,1,0},{0,0,5}}};
  tri = (triangle) {{{0,0,0},{3,0,0},{2,2,0}}};
  //This triangle lies in the same plane as a facet, however the triangle is not a facet. So not disjoint.
  if (tri_tet_disjoint(&tri, &tet))
    printf("ERROR5 IN DISJOINT TEST");

  tet = (tetra) {{{0,0,0},{2,0,0},{1,1,0},{0,0,5}}};
  tri = (triangle) {{{0,0,0},{2,0,0},{1,5,0}}};
  //This triangle lies in the same plane as a facet and shares an edge. Should not be disjoint
  if (tri_tet_disjoint(&tri, &tet))
    printf("ERROR5.5 IN DISJOINT TEST");

  tri = (triangle) {{{2,0,0},{4,0,0},{3,1,0}}};
  //Triangle lies in the same plane as a facet, shares an vertex. Should be disjoint
  if (!tri_tet_disjoint(&tri, &tet))
    printf("ERROR 5.55 IN DISJOINT TEST");

  tri = (triangle) {{{2,0,0},{3,0,0},{1,1,0}}};
  //Triangle lies in the same plane as a facet, shares an edge. Should be disjoint
  if (!tri_tet_disjoint(&tri, &tet))
    printf("ERROR 5.555 IN DISJOINT TEST");

  tet = (tetra) {{{0,0,0}, {2,0,0}, {1,1,0},{1,1,9}}};
  tri = (triangle) {{{0,2,0},{1,1,0},{2,2,0}}};
  //Triangle and tetra share a point, should be disjoint
  if (!tri_tet_disjoint(&tri, &tet))
    printf("ERROR. 5555 IN DISJOINT TEST");

  tet = (tetra) {{{0,0,0}, {1,0,0}, {0,1,0},{0,0,9}}};
  tri = (triangle) {{{1,0,0},{2,0,0},{2,1,0}}};
  //Triangle and tetra share a point, should be disjoint
  if (!tri_tet_disjoint(&tri, &tet))
    printf("ERROR. 5555 IN DISJOINT TEST");

  tetra t1 = (tetra) {{{5,5,5},{6,5,5},{5,6,5},{5,5,6}}};
  tetra t2 = (tetra) {{{5,5,5},{6,5,5},{5,6,5},{5,5,4}}};
  //These two tetrahedrons share a facet, and disjoint
  if (!tet_tet_disjoint(&t1,&t2))
    printf("ERROR 6 IN DISJOINT TEST");


  t1 = (tetra) {{{0,0,0},{16,0,0},{12,9,0},{10,5,9}}};
  t2 = (tetra) {{{0,0,0},{0,16,0},{12,9,0},{10,5,9}}};
  //These two tetrahedrons share a edge, base triangle lies in the same plane and are disjoint
  if (!tet_tet_disjoint(&t1,&t2))
    printf("ERROR 7 IN DISJOINT TEST");

  t1 = (tetra) {{{5, 5, 5}, {21, 5, 5},{17,14, 5},{15,10,14}}};
  t2 = (tetra) {{{ 5, 5, 5},{13,13,14},{ 5,15, 1},{17,14, 5}}};
  //These two tetrahedrons share an edge, base triangle lies in a different plane, and disjoint
  if (!tet_tet_disjoint(&t1,&t2))
    printf("ERROR 8 IN DISJOINT TEST");
  t1 = (tetra) {{{1,2,1},{0,2,1},{2,1,1},{0,1,0}}};
  t2 = (tetra) {{{1,0,0},{0,0,2},{2,1,0},{2,0,0}}};
  //These two tetrahedrons are disjoint, but not trivially
  if (!tet_tet_disjoint(&t1,&t2)) {
    printf("ERROR 9 IN DISJOINT TEST");
    return;
  }

  t1 = (tetra) {{{0,0,0}, {1,0,0}, {1,1,0},{1,0,1}}};
  t2 = (tetra) {{{0,0,0}, {1,0,0}, {1,1,0},{0,1,1}}};
  //These two tetrahedrons have a facet in the same plane, but not same facet, so not disjoint
  if (tet_tet_disjoint(&t1, &t2))
    printf("ERROR 10 IN DISJOINT TEST");
  //The same tetrahedrons.. Should intersect
  if (tet_tet_disjoint(&t1, &t2))
    printf("ERROR 11 IN DISJOINT TEST");
}

int main(void){
  tri_mem_list list;
  ptriangulation triang;
  list = mem_list_init(30, MEM_LIST_CUBE_SPARSE, MEM_LIST_FALSE);
  tetra tet = (tetra) {{{0,0,0},                                                              
			{14,0,0},                                                             
			{7,8,0},                                                              
			{9,5,7}}};
  triangle sides[4];
  tet_sides(&tet, sides);
  for (int i = 0; i < 4; i++)
    mem_list_cube_set(&list, sides + i);

  triang = triangulate_cube(&list, NULL, NULL);
  test_tetra_disjoint();
  /*
  test_mem_list_fund();
  
  test_prism(20);
  test_boundary(5);
  test_sym();
  test_tetra_normals();
  test_mem_list_tet();
  test_triangle_indices();
  test_tetra_disjunct();
  */
  test_tri_list();
}
