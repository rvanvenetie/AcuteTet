#pragma once
#include <string>
#include "triangle.h"
#include "domain.h"
#include "set.h"

/**
 *  Base class for filtering
 *  T should be of the 
 */
template <typename T>
class TriangleFilter {
  private:
    // reference to the triangle set
    T &_set;
    
    // reference to the domain
    const decltype(T::_domain) &_domain;

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
    inline bool conform(const Triangle<3> &triangle) const;

    // define sweep method
    bool sweep();

    // define the filter method
    bool filter();
};
