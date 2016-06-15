#pragma once
#include "domain.h"
class Square : public Domain<2, Square> {
  public:
    using Domain<2,Square>::_scale;

    Square(byte scale) : Domain<2,Square>(scale) {
      for (int x = 0; x < scale; x++)
        for (int y = 0; y < scale; y++) 
            _vertices.push_back({{x,y}});
      _size = _vertices.size();
    }

    inline vindex index(byte a,byte b) const { return a* _scale + b; }
    inline vindex index(const Vector<2> &idx) const { return index(idx[0], idx[1]); }

    // returns whether this edge lies on the boundary
    inline bool boundary(const Simplex<2,2> &edge) const 
    {
      return (edge.slice(0) == 0 || edge.slice(0) == _scale -1 ||
              edge.slice(1) == 0 || edge.slice(1) == _scale -1);
    }
};
