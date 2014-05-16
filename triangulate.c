#include "limits.h"
#include "triangle.h"
#include "tetraeder.h"
#include "triangulate.h"

void tetra_edges(ptetra tet, arr3 * edges) {
  //Edges of base triangle
  subArr3(tet->vertices[1], tet->vertices[0], edges[0]);
  subArr3(tet->vertices[2], tet->vertices[0], edges[1]);
  subArr3(tet->vertices[2], tet->vertices[1], edges[2]);
  //Now only need the three edges from base to apex 
  subArr3(tet->vertices[3], tet->vertices[0], edges[3]);
  subArr3(tet->vertices[3], tet->vertices[1], edges[4]); 
  subArr3(tet->vertices[3], tet->vertices[2], edges[5]);     
}
/*
 * Test for non-strict overlapping
 */
int tetra_disjoint_axis(ptetra t1, ptetra t2, arr3 axis) {
  int min_1 = INT_MAX, max_1 = INT_MIN;
  int min_2 = INT_MAX, max_2 = INT_MIN;
  for (int i = 0; i < 4; i++) {
    int dot1 = dotArr3(axis, t1->vertices[i]);
    if (dot1 < min_1) min_1 = dot1;
    if (dot1 > max_1) max_1 = dot1;

    int dot2 = dotArr3(axis, t2->vertices[i]);
    if (dot2 < min_2) min_2 = dot2;
    if (dot2 > max_2) max_2 = dot2;
  }
  //Now we know the scalar projections on the axis, test for overlap
  
  /*
   * Disjoint if range A is completely after range B, or that range A is completely before B.
   */
  return (min_1 >= max_2 || max_1 <= min_1);
}
/*
 * If disjoint, we return true, as this means we have found a separating plane
 */
int tetra_disjoint_axes(ptetra t1, ptetra t2, arr3 * axes, int len) {
  for (int i = 0; i < len; i++) {
    if (tetra_disjoint_axis(t1,t2,axes[i]))
      return 1;
  }
  return 0;
}
int tetra_tetra_disjoint(ptetra t1, ptetra t2) {
  arr3 normals[4];
  tetra_normals(t1, normals); //First check facets for separating axis
  if (tetra_disjoint_axes(t1,t2,normals, 4))
    return 1;
  tetra_normals(t2, normals);
  if (tetra_disjoint_axes(t1,t2,normals,4))
    return 1;
  
  //Indication how edge i is made up.. for example edge[0] = vertex[1] - vertex[0]
  int edges[6][2] = {{1,0},{2,0},{2,1}, //Base triangle
                   {3,0},{3,1},{3,2}}; //Edges from base to apex
  //Test all the edges
  arr3 edge_a, edge_b;
  for (int i = 0; i < 6; i++) {
    subArr3(t1->vertices[edges[i][0]], t1->vertices[edges[i][1]], edge_a);
    for (int j = 0; j < 6; j++)
    {
      subArr3(t2->vertices[edges[j][0]], t2->vertices[edges[j][1]], edge_b);
      //Find normal to plane parallel to edge_a and edge_b
      arr3 perp_axis;
      crossArr3(edge_a,edge_b, perp_axis);
      if (zeroArr3(perp_axis)) { //Edges lie in the same plane, try vector from a to b
        subArr3(t1->vertices[edges[i][0]], t2->vertices[edges[j][0]], edge_b);
        crossArr3(edge_a,edge_b, perp_axis);
        if (zeroArr3(perp_axis)) //Still zero, all points lie on the same line, not a separation axis, IS TRUE?
          continue;
      }     
      if (tetra_disjoint_axis(t1,t2, perp_axis))
        return 1;
    }
  }
  return 0; //No separation axis found, thus we are not disjoint
}

ptriangulation triangulate_cube_random(arr3 dim) {
  //ptriangulation result = calloc(sizeof(triangulation), 1);
  
}
