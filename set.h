#pragma once
#include <bitset>
#include <typeinfo>
#include "vector.h"
#include "triangle.h"
#include "domain.h"

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
    inline bool contains(const Triangle<3> &triang, const Vector<3> &apex) const
    {
      return (contains(triang[0], triang[1], apex) &&
              contains(triang[0], triang[2], apex) &&
              contains(triang[1], triang[2], apex));
    }

    // contains triangle determined by it's geometrical vertices
    inline bool contains(const Vector<3> &a, const Vector<3> &b, const Vector<3> &c) const
    {
      // determine indices
      tindex indices(super()->index(a,b,c));
      
      // sort indices
      sortindices(indices);

      // return whether we got the indices
      return contains(indices);
    }

    // delete
    inline void reset(const Triangle<3> &triangle) { super()->reset(triangle); }

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

class CubeTriangleSet : public TriangleSet<CubeTriangleSet, Cube> {
  private:
    // hold the size of the first dimension
    vindex _size = 0;
  public:
    // initalize cube triangle set
    CubeTriangleSet(byte scale, bool set=false) :
      TriangleSet<CubeTriangleSet, Cube>(scale),
      _size(scale * scale * scale)
    { 
      this->_name = "CubeTriangleSet";
      this->init(set); 
    }

    // define axis sizes for this data type
    inline vindex size() const { return _size; }
    inline vindex size(vindex a) const { return _size - a; }
    inline vindex size(vindex a, vindex b) const  { return _size - a - b; } 

    // define vertex to index
    inline vindex index(const Vector<3> &v) const { return _domain.index(v); }
    inline tindex index(const Vector<3> &a, const Vector<3> &b, const Vector<3> &c) const { return {{index(a), index(b), index(c)}}; }

    // reset triangle
    using TriangleSet<CubeTriangleSet, Cube>::reset;
    inline void reset(const Triangle<3> &triang) { reset(index(triang[0]), index(triang[1]), index(triang[2])); }

    // from file
    static CubeTriangleSet * fromFile(const std::string &filename);
};

// A memory set containing triangles for the cube, this filters
// out symmetries. So for each triangle, only one of its symmetry class is stored
class FundcubeTriangleSet : public TriangleSet<FundcubeTriangleSet, Cube> {
  private:
    // fundamental domain
    Fundcube _fund;

    // size cube
    vindex _cubesize;

    // size fund
    vindex _fundsize;

    // Fundcube vertices: fund idx -> cube pt
    vector<Vector<3>> _vertices;

    // index array: cube idx -> fund idx
    vector<vindex> _index;

    // symmetry index: cube idx -> sym idx
    vector<byte> _symindex;

    // converts : cube idx + sym -> fund_idx
    vector<vindex> _sympt;
  public:

    // initalize cube triangle set
    FundcubeTriangleSet(byte scale, bool set=false);

    // define axis sizes for this data type
    inline vindex size() const { return _fundsize; }
    inline vindex size(vindex a) const { return _cubesize - a; }
    inline vindex size(vindex a, vindex b) const  { return _cubesize - a - b; } 

    // define vertex to (fund) index
    inline vindex index(const Vector<3> &v) const { return _index[_domain.index(v)]; }
    // define triangle to (fund) indices
    inline tindex index(const Vector<3> &a, const Vector<3> &b, const Vector<3> &c) const 
    {
      size_t sym = _symindex[_domain.index(a)]; 
      return {{_sympt[_domain.index(a) + _cubesize*sym],
               _sympt[_domain.index(b) + _cubesize*sym],
               _sympt[_domain.index(c) + _cubesize*sym]}};
    }

    using TriangleSet<FundcubeTriangleSet, Cube>::reset;

    // delete all of its symmetries as well
    inline void reset(const Triangle<3> &triangle) 
    {
      for (size_t sym = 0; sym < 48; sym ++) 
        if (_sympt[_domain.index(triangle[0]) + _cubesize*sym] < _fundsize ||
            _sympt[_domain.index(triangle[1]) + _cubesize*sym] < _fundsize ||
            _sympt[_domain.index(triangle[2]) + _cubesize*sym] < _fundsize) 
        {
          tindex indices ={{_sympt[_domain.index(triangle[0]) + _cubesize*sym],
                            _sympt[_domain.index(triangle[1]) + _cubesize*sym],
                            _sympt[_domain.index(triangle[2]) + _cubesize*sym]}};
          sortindices(indices);
          reset(indices);
        }
    }

    // turns (fund) index to a vertex
    inline const Vector<3> &vertex(vindex a) const { return _vertices[a]; }

    // define the fund
    inline const Fundcube &fund() const { return _fund; }

    // from file
    static FundcubeTriangleSet * fromFile(const std::string &filename);
};
