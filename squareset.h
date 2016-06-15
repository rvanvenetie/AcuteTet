#pragma once
#include "set.h"
#include "square.h"

class SquareTriangleSet : public TriangleSet<SquareTriangleSet, Square> {
  private:
    // hold the size of the first dimension
    vindex _size = 0;

  public:
    // initalize cube triangle set
    SquareTriangleSet(byte scale, bool set=false) :
      TriangleSet<SquareTriangleSet, Square>(scale),
      _size(scale * scale)
    { 
      this->_name = "SquareTriangleSet";
      this->init(set); 
    }

    // define axis sizes for this data type
    inline vindex size() const { return _size; }
    inline vindex size(vindex a) const { return _size - a; }
    inline vindex size(vindex a, vindex b) const  { return _size - a - b; } 

    // define vertex to index
    inline vindex index(const Vector<2> &v) const { return _domain.index(v); }
    inline tindex index(const Vector<2> &a, const Vector<2> &b, const Vector<2> &c) const { return {{index(a), index(b), index(c)}}; }

    // reset triangle
    using TriangleSet<SquareTriangleSet, Square>::reset;
    inline void reset(const Triangle<2> &triang) { reset(index(triang[0]), index(triang[1]), index(triang[2])); }

    // from file
    static SquareTriangleSet * fromFile(const std::string &filename);
};
