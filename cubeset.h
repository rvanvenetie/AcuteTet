#pragma once
#include "set.h"
#include "cube.h"
#include "squareset.h"
// S is the implementor of the set
template <typename S>
class CubeTSet : public S
{
  private:
    // Cube domain
    Cube _domain;
  public:

    // Constructor from file
    CubeTSet(std::string filename);

    // Initalize cube triangle set
    CubeTSet(byte scale, bool set=false) : _domain(scale)
    { 
      // axis length
      this->init({_domain.size(), _domain.size(), _domain.size()}, scale, set);
      this->_name = "CubeTSet";
    }

    // define vertex to index
    inline vindex index(const Vector<3> &v) const { return _domain.index(v); }
    inline tindex index(const Vector<3> &a, const Vector<3> &b, const Vector<3> &c) const { return {{index(a), index(b), index(c)}}; }

    // turns index to a vertex
    inline const Vector<3> &vertex(vindex a) const { return _domain[a]; }

    // reset triangle
    inline void reset(const Triangle<3> &triang) { reset(index(triang[0]), index(triang[1]), index(triang[2])); }
};

// A memory set containing triangles for the cube, this filters
// out symmetries. So for each triangle, only one of its symmetry class is stored
template <typename S>
class FCubeTSet : public S {
  private:
    // Cube domain
    Cube _domain;

    // fundamental domain
    Fundcube _fund;

    // sizes domains
    vindex _cubesize;
    vindex _fundsize;

    // FCube vertices: fund idx -> cube pt
    vector<Vector<3>> _vertices;

    // index array: cube idx -> fund idx
    vector<vindex> _index;

    // symmetry index: cube idx -> sym idx
    vector<byte> _symindex;

    // converts : cube idx + sym -> fund_idx
    vector<vindex> _sympt;

    // inits all the above variables to the correct values
    void init(byte scale);

    // legacy load
    void legacyload(const std::string &filename);
  public:
    const static int dim = 3;

    // from file
    FCubeTSet(const std::string &filename, bool legacy=false);

    // initalize fund cube triangle set
    FCubeTSet(byte scale, bool set=false);

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

    using S::reset;
    using S::sortindices;
    using S::contains;
    using S::set;


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
      tindex indices(index(a,b,c));
      
      // sort indices
      sortindices(indices);

      // return whether we got the indices
      return contains(indices);
    }

    // contains triangl
    inline bool contains(const Triangle<3> tri) const { return contains(tri[0], tri[1], tri[2]); }

    // turns (fund) index to a vertex
    inline const Vector<3> &vertex(vindex a) const { return _vertices[a]; }

    // define the fund
    inline const Fundcube &fund() const { return _fund; }

    // define the fund
    inline const Cube &domain() const { return _domain; }
};
