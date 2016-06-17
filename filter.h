#pragma once
#include <string>
#include "triangle.h"
#include "domain.h"
#include "cube.h"
#include "set.h"
#include "cubeset.h"

/**
 *  Base class for filtering
 *  T should be the class used to store the triangles
 */
template <typename T>
class TriangleFilter {
  private:
    // reference to the triangle set
    T &_set;
    
    // reference to the domain
    const decltype(((T*) nullptr)->domain()) &_domain;

    // store file name
    std::string _filename;
    std::string _tmpfile;

    // store interval
    size_t _interval;
  public:
    
    TriangleFilter(T &set, std::string filename = "", std::string tmpfile = "", size_t interval = 0) :
      _set(set),
      _domain(set.domain()),
      _filename(filename),
      _tmpfile(tmpfile),
      _interval(interval) {}

    // define conformity checker for triangle given by vertices (a,b,c)
    inline bool valid(const Triangle<T::dim> &triangle) const;

    // define sweep method
    bool sweep();

    // define the filter method
    bool filter();

    // boundary facets, works for all 3D-dimensional sets
    //   intersects with a plane perpendicular to axis at height point, and returns all such points
    SquareTSet boundaryfacets(byte axis, byte point) const;
    // intersects with one of the 6 sides
    SquareTSet boundaryfacets(byte side) const {
      int axis[] = {0,1,2,0,1,2};
      int pt[]   = {0,0,0, _set.scale()-1,_set.scale()-1, _set.scale()-1};
      return boundaryfacets(axis[side], pt[side]);
    }

    SquareTSet boundaryfacets() const { return boundaryfacets(0); }
};
