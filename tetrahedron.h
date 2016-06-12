#pragma once
#include "vector.h"
#include "triangle.h"
#include "polygon.h"

struct Tetrahedron  : public Polygon<3,4> 
{
  using Polygon<3,4>::_vertices;


  // move constructor
  inline Tetrahedron(Polygon<3,4> &&other)  : Polygon<3,4>(other) {}

  // constructor
  inline Tetrahedron(const Triangle<3> &triangle, const Vector<3> & pt) : Polygon<3,4>({triangle[0], triangle[1], triangle[2], pt}) {}

  // test whether this tetrahedron is acute
  inline static bool acute(const Tetrahedron &tet) 
  {
    Vector<3> edges[5];
    Vector<3> normals[4];

    edges[0] = tet[1] - tet[0]; // b - a
    edges[1] = tet[2] - tet[0]; // c - a
    edges[2] = tet[3] - tet[0]; // d - a
    edges[3] = tet[2] - tet[1]; // c - b
    edges[4] = tet[3] - tet[1]; // d - b

    normals[0] = cross(edges[4] , edges[3]);
    normals[1] = cross(edges[1] , edges[2]);
    normals[2] = cross(edges[2] , edges[0]);
    normals[3] = cross(edges[0] , edges[1]);

    return (dot(normals[0], normals[1]) < 0 &&
            dot(normals[0], normals[2]) < 0 &&
            dot(normals[0], normals[3]) < 0 &&
            dot(normals[1], normals[2]) < 0 &&
            dot(normals[1], normals[3]) < 0 &&
            dot(normals[2], normals[3]) < 0);
  }
  inline bool acute() const { return acute(*this); }

  // optimized version, test whether tetrahedron is acute, given by triangle and apex
  
  inline static bool acute(const Triangle<3> &triangle, const Vector<3> &pt)
  {
    Vector<3> edges[5];
    Vector<3> normals[4];

    // three edges of the tetra
    edges[0] = triangle[0] - pt;
    edges[1] = triangle[1] - pt;
    edges[2] = triangle[2] - pt;

    // all facets must be acute, cheap test to rule this tetra out right now
    if (dot(edges[0], edges[1]) <= 0 || dot(edges[0], edges[2]) <= 0 || dot(edges[1],edges[2]) <= 0)
      return false;

    normals[2] = cross(edges[2], edges[0]);
    normals[3] = cross(edges[0], edges[1]);

    if (dot(normals[2], normals[3]) >= 0)
      return false;

    normals[1] = cross(edges[1], edges[2]);
    if (dot(normals[1], normals[2]) >= 0 || dot(normals[1], normals[3]) >= 0)
      return false;

    edges[3] = triangle[1] - triangle[0];
    edges[4] = triangle[2] - triangle[0];

    normals[0] = cross(edges[4], edges[3]);

    return (dot(normals[0], normals[1]) < 0 &&
            dot(normals[0], normals[2]) < 0 &&
            dot(normals[0], normals[3]) < 0);
  }

};
