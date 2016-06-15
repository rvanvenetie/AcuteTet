#pragma once
#include "domain.h"

class Cube : public Domain<3, Cube> {
  public:
    using Domain<3,Cube>::_scale;

    Cube(byte scale) : Domain<3,Cube>(scale) {
      for (int x = 0; x < scale; x++)
        for (int y = 0; y < scale; y++) 
          for (int z = 0; z < scale; z++) {
            _vertices.push_back({{x,y,z}});
          }
      _size = _vertices.size();
    }

    inline vindex index(byte a,byte b,byte c) const { return a*_scale * _scale + b*_scale + c; }
    inline vindex index(const Vector<3> &idx) const { return index(idx[0], idx[1], idx[2]); }

    // returns whether this triangle lies on the boundary
    inline bool boundary(const Simplex<3,3> &triang) const 
    {
      return (triang.slice(0) == 0 || triang.slice(0) == _scale -1 ||
              triang.slice(1) == 0 || triang.slice(1) == _scale -1 ||
              triang.slice(2) == 0 || triang.slice(2) == _scale - 1);
    }
};

class Fundcube : public Domain<3,Fundcube> {
  public:
    // vertices
    using Domain<3,Fundcube>::_vertices;
    
    Fundcube(byte scale)  : Domain<3,Fundcube>(scale) {
      int M = (scale-1) / 2; //has points 0..scale-1, automatically floored. As point 3.5 is inside
      for (int z = 0; z <= M; z++)
        //We now have triangle (z,z,z) - (z,m,z) - (m,m,z)
        for (int x = z; x <= M; x++) //Loop over x-axis
          for (int y = z; y <= x; y++) {
            _vertices.push_back({{x,y,z}});
          }
      _size = _vertices.size();
    }

    // does the domain contain this vertex?
    inline bool contains(const Vector<3> &vertex) const
    {
      int M = (_scale-1) / 2; //has points 0..scale-1, automatically floored. As point 3.5 is inside
      return  0 <= vertex[2]         && vertex[2] <= M &&
              vertex[2] <= vertex[0] && vertex[0] <= M &&
              vertex[2] <= vertex[1] && vertex[1] <= vertex[0];
    }


    // returns the symmetry index for a given vertex
    inline byte symindex(const Vector<3> &vertex)  const
    {
      for (byte sym = 0; sym < 48; sym ++)
        if (contains(symmetry(vertex, sym))) return sym;
      return 0;
    }

    // applies symmetry on the given vertex
    inline Vector<3> symmetry(const Vector<3> &vertex, byte sym) const 
    {
      #define mirror {t = x;x = y;y = t;}    
      #define face0 {}
      #define face1 {t = y;y = _scale - 1 - z;z = t;}
      #define face2 {t = x;x = _scale - 1 - z;z = t;}
      #define face3 {y = _scale - 1 - y;z = _scale - 1 - z;}
      #define face4 {t = y;y = z;z = _scale - 1 - t;}
      #define face5 {t = x;x = z;z = _scale - 1 - t;}
      #define rot0 {}
      #define rot90 {t = x; x = y; y = _scale - 1 - t;}
      #define rot180 {x = _scale - 1 - x;y = _scale -1 - y;}
      #define rot270 {t = x;x = _scale - 1 - y;y = t;}
      byte t, x = vertex[0], y = vertex[1], z = vertex[2];
      switch (sym) {
        case 0:  face0; rot0;   break;
        case 1:  face0; rot90;  break;
        case 2:  face0; rot180; break;
        case 3:  face0; rot270; break;
        case 4:  face1; rot0;   break;
        case 5:  face1; rot90;  break;
        case 6:  face1; rot180; break;
        case 7:  face1; rot270; break;
        case 8:  face2; rot0;   break;
        case 9:  face2; rot90;  break;
        case 10: face2; rot180; break;
        case 11: face2; rot270; break;
        case 12: face3; rot0;   break;
        case 13: face3; rot90;  break;
        case 14: face3; rot180; break;
        case 15: face3; rot270; break;
        case 16: face4; rot0;   break;
        case 17: face4; rot90;  break;
        case 18: face4; rot180; break;
        case 19: face4; rot270; break;
        case 20: face5; rot0;   break;
        case 21: face5; rot90;  break;
        case 22: face5; rot180; break;
        case 23: face5; rot270; break;
        case 24: mirror; face0; rot0;   break;
        case 25: mirror; face0; rot90;  break;
        case 26: mirror; face0; rot180; break;
        case 27: mirror; face0; rot270; break;
        case 28: mirror; face1; rot0;   break;
        case 29: mirror; face1; rot90;  break;
        case 30: mirror; face1; rot180; break;
        case 31: mirror; face1; rot270; break;
        case 32: mirror; face2; rot0;   break;
        case 33: mirror; face2; rot90;  break;
        case 34: mirror; face2; rot180; break;
        case 35: mirror; face2; rot270; break;
        case 36: mirror; face3; rot0;   break;
        case 37: mirror; face3; rot90;  break;
        case 38: mirror; face3; rot180; break;
        case 39: mirror; face3; rot270; break;
        case 40: mirror; face4; rot0;   break;
        case 41: mirror; face4; rot90;  break;
        case 42: mirror; face4; rot180; break;
        case 43: mirror; face4; rot270; break;
        case 44: mirror; face5; rot0;   break;
        case 45: mirror; face5; rot90;  break;
        case 46: mirror; face5; rot180; break;
        case 47: mirror; face5; rot270; break;         
      }
    return {{x,y,z}};
    }
};
