#pragma once
#include <bitset>
#include <typeinfo>
#include "vector.h"
#include "triangle.h"
#include "domain.h"
#include "cube.h"
#include "square.h"

// A memory set containing all the triangles
// -  T is the type of the implementation class
// -  D is the type of the domain
template <typename T, typename D>
class TriangleSet {
  private:
    // short to parent
    inline const T * super() const { return static_cast<const T*>(this); }
    inline T * super() { return static_cast<T*>(this); }

  protected:
    // hidden constructor
    TriangleSet(byte scale) : _scale(scale), _domain(scale) {}

    // initalization method
    void init(bool set);


  public:
    //boost::dynamic_bitset<>  ** _data2 = nullptr;
    byte  *** _data = nullptr;         //Actual data array
    byte _scale;                       //Scale
    D _domain;                         //Domain
    const char *_name = "TriangleSet"; //Class name
    static const int dim = D::dim;
    

    // constructor
    TriangleSet(byte scale, bool set) : _scale(scale), _domain(scale) { init(set); }

    // destructor
    ~TriangleSet();

    // remove copy constructor
    TriangleSet(const TriangleSet &rhs) = delete;
    TriangleSet(TriangleSet &rhs) = delete;

    // return size of the axis'
    inline vindex size() const { return super()->size(); }
    inline vindex size(vindex a) const { return super()->size(a); }
    inline vindex size(vindex a, vindex b) const { return super()->size(a,b); } 

    // getters for scale and domain
    inline byte scale() const { return _scale; }
    inline const D & domain() const { return _domain; }

    // contains, set, reset for vertices (a,b,c)
    // 
    inline bool contains(vindex a, vindex b, vindex c) const { return _data[a][b][c/ 8] & (1 << (c % 8)); }
    inline void set(vindex a, vindex b, vindex c) { _data[a][b][c/8] |= 1 << (c % 8); }
    inline void reset(vindex a, vindex b, vindex c) { _data[a][b][c/8] &= ~(1 << (c%8)); }
    inline void reset(const tindex &tri) { reset(tri[0], tri[1], tri[2]); }
    /*
    inline bool contains(vindex a, vindex b, vindex c) const { return _data2[a][b].test(c); }
    inline void set(vindex a, vindex b, vindex c) { _data2[a][b].set(c); } 
    inline void reset(vindex a, vindex b, vindex c) { _data2[a][b].reset(c); }
    */

    // contains triangle by indeices
    inline bool contains(const tindex &v) const { return contains(v[0], v[1], v[2]); }
    inline bool operator()(vindex a, vindex b, vindex c) const { return contains(a,b,c); }
    inline void sortindices( tindex &indices) const 
    {
      // sort
      if (indices[0] > indices[1]) swap(indices[0],indices[1]);
      if (indices[1] > indices[2]) swap(indices[1],indices[2]);
      if (indices[0] > indices[1]) swap(indices[0],indices[1]);

      // subtract
      indices[2] = indices[2] - indices[1];
      indices[1] = indices[1] - indices[0];  
    }

    // returns whether it contains all the facets of tetrahedrons given by triangle and apex. Does not check triangle itself
    inline bool contains(const Triangle<D::dim> &triang, const Vector<D::dim> &apex) const
    {
      return (contains(triang[0], triang[1], apex) &&
              contains(triang[0], triang[2], apex) &&
              contains(triang[1], triang[2], apex));
    }

    // contains triangle determined by it's geometrical vertices
    inline bool contains(const Vector<D::dim> &a, const Vector<D::dim> &b, const Vector<D::dim> &c) const
    {
      // determine indices
      tindex indices(super()->index(a,b,c));
      
      // sort indices
      sortindices(indices);

      // return whether we got the indices
      return contains(indices);
    }

    // delete
    inline void reset(const Triangle<D::dim> &triangle) { super()->reset(triangle); }

    // returns the number of triangles starting with index i,j
    size_t count(vindex i, vindex j) const;

    // returns the amount of set triangles
    size_t count() const;

    // returns the memory used
    size_t memory() const;
    
    // outputs debug information, returns the number of triangles
    size_t print() const;

    // file options
    enum FileOptions { FULL, SPARSE };

    // to file
    bool toFile(const std::string &filename, FileOptions options = FULL) const;

    // from filestream
    bool fromFile(ifstream &file);
};
