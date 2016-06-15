#pragma once
#include "set.h"
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
