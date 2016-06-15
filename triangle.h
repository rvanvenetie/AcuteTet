#pragma once
#include "vector.h"
#include "simplex.h"

struct tindex 
{
  vindex idx[3];
  vindex const &operator[] (byte i) const { return idx[i]; }
  vindex &operator[] (byte i) { return idx[i]; }
  void print() const { cout << "(" << idx[0] << "," << idx[1] << "," << idx[2] << ")" << endl; }
};

template <byte _dim>
struct Triangle : public Simplex<_dim,3> 
{
  // move constructors
  inline Triangle(Simplex<_dim,3> &&other)  : Simplex<_dim,3>(move(other)) {}

  inline Simplex<_dim,3> edges() const;
  inline bool static acute(const Simplex<_dim,3> &edges);
  inline bool acute() const { return acute(edges()); }
};

template<>
inline Simplex<3,3> Triangle<3>::edges() const
{
  return {{_vertices[1] - _vertices[0], _vertices[2] - _vertices[0], _vertices[2] - _vertices[1]}};
}

template<>
inline Simplex<2,3> Triangle<2>::edges() const
{
  return {{_vertices[1] - _vertices[0], _vertices[2] - _vertices[0], _vertices[2] - _vertices[1]}};
}

template <byte _dim>
inline bool Triangle<_dim>::acute(const Simplex<_dim,3> &edges)
{
  return (dot(edges[0], edges[1]) > 0 && dot(edges[1], edges[2]) > 0 && dot(edges[0], edges[2]) < 0);
}

