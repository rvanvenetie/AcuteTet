#pragma once
#include "vector.h"

struct Tetrahedron {
  Vector<3> vertices[4];

  // direct access to the vertices
  Vector<3> & operator[](byte idx) { return vertices[idx]; }

  // print the tetrahedron
  void print() {
	printf("[[%d,%d,%d],\n",vertices[0][0],vertices[0][1],vertices[0][2]);
	printf(" [%d,%d,%d],\n",vertices[1][0],vertices[1][1],vertices[1][2]);
	printf(" [%d,%d,%d],\n",vertices[2][0],vertices[2][1],vertices[2][2]);
	printf(" [%d,%d,%d]]\n",vertices[3][0],vertices[3][1],vertices[3][2]);
  }

  static Tetrahedron random(byte dim) {
    return {Vector<3>::random(dim),
            Vector<3>::random(dim),
            Vector<3>::random(dim),
            Vector<3>::random(dim)};
  }
  inline int acute() const {
    Vector<3> edges[5];
    Vector<3> normals[4];

    edges[0] = vertices[1] - vertices[0]; // b - a
    edges[1] = vertices[2] - vertices[0]; // c - a
    edges[2] = vertices[3] - vertices[0]; // d - a
    edges[3] = vertices[2] - vertices[1]; // c - b
    edges[4] = vertices[3] - vertices[1]; // d - b

    normals[0] = edges[4] * edges[3];
    normals[1] = edges[1] * edges[2];
    normals[2] = edges[2] * edges[0];
    normals[3] = edges[0] * edges[1];

    return (dot(normals[0], normals[1]) < 0 &&
            dot(normals[0], normals[2]) < 0 &&
            dot(normals[0], normals[3]) < 0 &&
            dot(normals[1], normals[2]) < 0 &&
            dot(normals[1], normals[3]) < 0 &&
            dot(normals[2], normals[3]) < 0);
  }
};
