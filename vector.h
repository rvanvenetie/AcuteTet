#pragma once
#include <iostream>

using namespace std;
typedef unsigned char byte;
template <byte len>
struct Vector{
    int entries[len];

    // direct access to the entries
    int & operator[](byte idx) { return entries[idx]; }

    // random vector
    static Vector<len> random(byte dim) {
      Vector<len> result;
      for (byte i =0; i < len; i++)
        result.entries[i] = rand() % dim;

      return result;
    }

    // print vector
    void print() const {
      for (byte i =0; i < len; i++) 
        cout << ((i == 0)? "[": ",") << entries[i];
      cout << "]\n";
    }

    // zero vector?
    inline bool zero() const {
      for (byte i =0; i < len; i++) 
        if (entries[i]) return false;
      return true;
    }

    // do two objects equal?
    inline bool operator==(const Vector<len>& rhs) const {
      for (byte i =0; i < len; i++)
        if (entries[i] != rhs.entries[i]) return false;
      return true;
    }

    // subtract two vectors
    inline Vector<len> operator -(const Vector<len> &rhs) const {
      Vector<len> result;
      for (byte i = 0; i < len; i++)
        result.entries[i] = entries[i] - rhs.entries[i];
      return result;
    }

    // dot product
    inline int dot(const Vector<len> &rhs) const {
      int result = 0;
      for (byte i =0; i < len; i++)
        result += entries[i] * rhs.entries[i];
      return result;
    }

};

// add cross product for three lenensional operators
inline Vector<3> operator *(const Vector<3>& lhs, const Vector<3> rhs){ 
  return {lhs.entries[1] * rhs.entries[2] - lhs.entries[2] * rhs.entries[1],
          lhs.entries[2] * rhs.entries[0] - lhs.entries[0] * rhs.entries[2],
          lhs.entries[0] * rhs.entries[1] - lhs.entries[1] * rhs.entries[0]};
}
// add dot product to the name space
template <byte len>
inline int dot (const Vector<len> &lhs, const Vector<len> &rhs) {
      int result = 0;
      for (byte i =0; i < len; i++)
        result += lhs.entries[i] * rhs.entries[i];
      return result;
    }

