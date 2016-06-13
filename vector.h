#pragma once
#include <iostream>

using namespace std;

typedef unsigned char byte;
typedef unsigned short vindex;

template <byte _dim>
struct Vector{
    int entries[_dim];

    // direct access to the entries
    inline int & operator[](byte idx) { return entries[idx]; }
    inline int const & operator[](byte idx) const { return entries[idx]; }

    // random vector with points  in the mesh grid [0,scale]^_dim
    static Vector<_dim> random(byte scale) {
      Vector<_dim> result;
      for (byte i =0; i < _dim; i++)
        result.entries[i] = rand() % scale;

      return result;
    }

    // print vector
    void print() const {
      for (byte i =0; i < _dim; i++) 
        cout << ((i == 0)? "[": ",") << entries[i];
      cout << "]\n";
    }

    // zero vector?
    inline bool zero() const {
      for (byte i =0; i < _dim; i++) 
        if (entries[i]) return false;
      return true;
    }

    // do two objects equal?
    inline bool operator==(const Vector<_dim>& rhs) const {
      for (byte i =0; i < _dim; i++)
        if (entries[i] != rhs.entries[i]) return false;
      return true;
    }

    // does every index equal?
    inline bool operator==(byte rhs) const 
    {
      for (byte i =0; i < _dim; i++)
        if (entries[i] != rhs) return false;
      return true;
    }

    // subtract two vectors
    /*
    inline Vector<_dim> operator -(const Vector<_dim> &rhs) const {
      Vector<_dim> result;
      for (byte i = 0; i < _dim; i++)
        result.entries[i] = entries[i] - rhs.entries[i];
      return result;
    }
    */

    // dot product
    inline int dot(const Vector<_dim> &rhs) const {
      int result = 0;
      for (byte i =0; i < _dim; i++)
        result += entries[i] * rhs.entries[i];
      return result;
    }

};

// add specific three dimensional operators; for speed.

// subtract two vectors for three dimensional operators
inline Vector<3> operator -(const Vector<3>& lhs, const Vector<3> &rhs){ 
  return {{lhs.entries[0] -rhs.entries[0], lhs.entries[1] - rhs.entries[1], lhs.entries[2] - rhs.entries[2]}};
}

// add cross product for three dimensional vectors
inline Vector<3> cross(const Vector<3>& lhs, const Vector<3> &rhs) {
  return {{lhs.entries[1] * rhs.entries[2] - lhs.entries[2] * rhs.entries[1],
          lhs.entries[2] * rhs.entries[0] - lhs.entries[0] * rhs.entries[2],
          lhs.entries[0] * rhs.entries[1] - lhs.entries[1] * rhs.entries[0]}};
}
inline void cross(const Vector<3>& lhs, const Vector<3> &rhs, Vector<3> &result) {
  result = {{lhs.entries[1] * rhs.entries[2] - lhs.entries[2] * rhs.entries[1],
          lhs.entries[2] * rhs.entries[0] - lhs.entries[0] * rhs.entries[2],
          lhs.entries[0] * rhs.entries[1] - lhs.entries[1] * rhs.entries[0]}};
}

// dot product for two dimensional vectors
inline int dot(const Vector<3> &lhs, const Vector<3> &rhs) { return lhs.entries[0] * rhs.entries[0] + lhs.entries[1] * rhs.entries[1] + lhs.entries[2] * rhs.entries[2]; }
