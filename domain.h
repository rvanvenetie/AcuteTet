#pragma once
#include <vector>
#include "vector.h"
#include "triangle.h"

using namespace std;

// Domain is a class that holds points.
//   Scale : Indicates the finess of the grid
//   Dim   : Indicates the dimension of the points
template <byte _dim, typename T>
class Domain {
  private:
    // short to parent
    inline const T * super() const { return static_cast<const T*>(this); }

  public:
    // static int exposing the dimension
    static const byte dim = _dim;
    // vector of the vertices
    vector<Vector<_dim>> _vertices;
    // scale of the grid
    byte _scale;
    // amount of points
    size_t _size;

    // Constructor
    Domain(byte scale) : _scale(scale) {} 

    // direct access to the entries
    inline const Vector<_dim> &operator[](vindex idx) const { return _vertices[idx]; }

    // member variables getters
    inline const  size_t &size() const { return _size; }
    inline size_t scale() const { return _scale; }

    // facet on the boundary
    inline bool boundary(const Simplex<_dim,_dim> &facet) const { return super()->boundary(facet); }
};
