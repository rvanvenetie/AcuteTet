#pragma once
#include "square.h"
#include "set.h"

class SquareTSet : public TFullSet<2> 
{
  private:
    // alias
    using S= TFullSet<2>;
    // Domain
    Square _domain;

  public:
    // initalize cube triangle set
    SquareTSet(byte scale, bool set=false) : _domain(scale)
    { 
      // axis length
      S::init({(vindex) _domain.size(),(vindex) _domain.size(),(vindex) _domain.size()}, scale, set);
      S::_name = "SquareTSet";
    }

    // Expose base  definitions, i.e. implementations using the index
    using S::contains;
    using S::reset;
    using S::set;

    // define vertex to index
    inline vindex index(const Vector<2> &v) const { return _domain.index(v); }
    inline tindex index(const Vector<2> &a, const Vector<2> &b, const Vector<2> &c) const { return {{index(a), index(b), index(c)}}; }

    // reset triangle
    inline void reset(const Triangle<2> &triang) { reset(index(triang[0]), index(triang[1]), index(triang[2])); }

    //does the set contain the triangle?
    inline bool contains(const Vector<2> &a, const Vector<2> &b, const Vector<2> &c) const
    {
      // determine indices
      tindex indices(index(a,b,c));
      
      // sort indices
      sortindices(indices);

      // return whether we got the indices
      return S::contains(indices);
    }

    // return domain
    const Square &domain() const { return _domain; }
};
