#pragma once
#include "vector.h"


template <byte _dim, byte _n>
struct Polygon {
  Vector<_dim> _vertices[_n];

  // direct access to the vertices
  Vector<_dim> & operator[](byte idx) { return _vertices[idx]; }
  Vector<_dim> const & operator[](byte idx) const { return _vertices[idx]; }

  // return a vector containing the coordinates of the specific axis
  Vector<_n> slice(byte axis) const { 
    Vector<_n> result;
    for (size_t i = 0; i < _n; i++) result[i] = _vertices[i][axis];
    return result;
  }

  // print the polygon
  void print() const {
    for (size_t i = 0; i < _n; i++) {
      cout << ((i == 0)?"[[":" [");
      for (size_t j = 0; j  < _dim; j++) {
        if (j !=0) cout << ",";
        cout << _vertices[i][j];
      }
      cout << ((i == _n-1)?"]]" : "],") << endl;
    }
  }

  static Polygon<_dim, _n> random(byte scale) 
  {
    Polygon<_dim, _n> result;
    for (size_t i = 0; i < _n; i++) result[i] = Vector<_dim>::random(scale);
    return result;
  }
};
