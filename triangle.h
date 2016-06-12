#pragma once
#include "vector.h"
#include "polygon.h"


struct tindex 
{
  vindex idx[3];
  vindex const &operator[] (byte i) const { return idx[i]; }
  vindex &operator[] (byte i) { return idx[i]; }
  void print() const { cout << "(" << idx[0] << "," << idx[1] << "," << idx[2] << ")" << endl; }
};

template <byte _dim>
struct Triangle : public Polygon<_dim,3> 
{
  // move constructors
  inline Triangle(Polygon<_dim,3> &&other)  : Polygon<_dim,3>(other) {}
  inline Triangle(Vector<_dim> &&a, Vector<_dim> &&b, Vector<_dim> &&c) : Polygon<_dim,3>({{a,b,c}}) {}
  // copy constructor
  inline Triangle(Vector<_dim> const &a, Vector<_dim> const &b, Vector<_dim> const &c) : Polygon<_dim,3>({{a,b,c}}) {}

  inline Polygon<_dim,3> edges() const;
  inline bool acute(const Polygon<_dim,3> &edges) const;
  inline bool acute() const 
  {
    return acute(edges());
  }
};

template<>
inline Polygon<3,3> Triangle<3>::edges() const
{
  return {{_vertices[1] - _vertices[0], _vertices[2] - _vertices[0], _vertices[2] - _vertices[1]}};
}

template <>
inline bool Triangle<3>::acute(const Polygon<3,3> &edges) const
{
  return (dot(edges[0], edges[1]) > 0 && dot(edges[1], edges[2]) > 0 && dot(edges[0], edges[2]) < 0);
}

