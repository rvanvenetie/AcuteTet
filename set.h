#pragma once
#include <bitset>
#include <typeinfo>
#include "vector.h"
#include "triangle.h"
#include "domain.h"

// Abstract class for Triangle Set. Uses CRTP.
// - T is the implementor 
// - D is the dimension
template <typename T, byte D>
class TSet 
{
  private:
    // short to implementationparent
    inline const T * super() const { return static_cast<const T*>(this); }
    inline T * super() { return static_cast<T*>(this); }

  public:
    const char * _name = "TSet";
    const static byte dim =  D;
    byte _scale;                 //scale used

    TSet() {}

    // implicit move constructor
    TSet(TSet &&) = default;

    // remove copy constructor
    TSet(const TSet &rhs) = delete;
    TSet(TSet &rhs) = delete;

    // getters
    byte scale() const { return _scale; }
    const char * name() const { return _name; }

    // main functionality, to be implemented by parent
    inline bool contains(vindex a, vindex b, vindex c) const { return super()->contains(a,b,c); }
    inline void set(vindex a, vindex b, vindex c) { return super()->set(a,b,c); }
    inline void reset(vindex a, vindex b, vindex c) { return super()->reset(a,b,c); }

    // return size of the axis, to be implemeted by parent
    inline vindex size() const { return super()->size(); }
    inline vindex size(vindex a) const { return super()->size(a); }
    inline vindex size(vindex a, vindex b) const { return super()->size(a,b); } 

    // some additional functions 
    inline void reset(const tindex &tri) { reset(tri[0], tri[1], tri[2]); }
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

    // returns the number of triangles starting with index i,j
    size_t count(vindex i, vindex j) const { return super()->count(i,j); }

    // returns whether this axis is empty
    bool empty(vindex i, vindex j) const { return super()->empty(i,j); }

    // returns the amount of set triangles
    size_t count() const;

    // returns the memory used
    size_t memory() const;
    
    // outputs debug information, returns the number of triangles
    size_t print() const;

    // to file
    bool toFile(const std::string &filename, bool sparse = false) const { return super()->toFile(filename, sparse); }
};


// Implementor of the abstract set class.
// Here we use FULL memory.
// - D is the dimension for which we hold triangles
template <byte D>
class TFullSet : public TSet<TFullSet<D>, D> {
  private:
    tindex _size;                    //Axis sizes

    using S=TSet<TFullSet<D>,D>;     //Alias to super
  protected:
    // Hidden constructor
    TFullSet() {}

    // Hidden initalization method
    void init(tindex size, byte scale, bool set);
    // constructor from file
    void init(ifstream &file, tindex size, byte scale);
  public:
    using S::_name;
    using S::_scale;
    byte  *** _data = nullptr;         //Actual data array

    // implicit move constructor
    TFullSet(TFullSet &&) = default;

    // constructor
    TFullSet(tindex size, byte scale, bool set) { init(size, scale, set); }

    // constructor from file
    TFullSet(ifstream &file, tindex size, byte scale) { init(file, size, scale); }

    // deconstructor
    ~TFullSet();

    // helper functions from parent
    using S::contains;
    using S::set;
    using S::reset;

    /**
     *  Necessary implementations
     */
    inline bool contains(vindex a, vindex b, vindex c) const { return _data[a][b][c/ 8] & (1 << (c % 8)); }
    inline void set(vindex a, vindex b, vindex c) { _data[a][b][c/8] |= 1 << (c % 8); }
    inline void reset(vindex a, vindex b, vindex c) { _data[a][b][c/8] &= ~(1 << (c%8)); }

    inline vindex size() const { return _size[0]; } 
    inline vindex size(vindex a) const { return _size[1] - a; }
    inline vindex size(vindex a, vindex b) const { return _size[2] - a - b; } 

    size_t count(vindex i, vindex j) const;
    bool empty(vindex i, vindex j) const;
    bool toFile(const std::string &filename, bool sparse = false) const;
};

// Implementor of the abstract set class.
// Here we use SPARSE memory.
// - D is the dimension for which we hold triangles
template <byte D>
class TSparseSet : public TSet<TSparseSet<D>, D> {
  private:
    tindex _size;                    //Axis sizes

    using S=TSet<TSparseSet<D>,D>;     //Alias to super
  protected:
    // hidden constructor
    TSparseSet() {};

    // Hidden initalization method
    void init(tindex size, byte scale, bool set);
    // constructor from file
    void init(ifstream &file, tindex size, byte scale);

  public:
    using S::_name;
    using S::_scale;
    byte  *** _data = nullptr;         //Actual data array

    // implicit move constructor
    TSparseSet(TSparseSet &&) = default;

    // constructor
    TSparseSet(tindex size, byte scale, bool set) { init(size, scale, set); }

    // constructor from file
    TSparseSet(ifstream &file, tindex size, byte scale) { init(file, size, scale); }

    // deconstructor
    ~TSparseSet();

    // checks whether the axies [a][b] exists
    inline bool exist(vindex a, vindex b) const { return _data[a][b] != nullptr; }

    // creates the axis [a][b] if needed
    inline void create(vindex a, vindex b) { 
      if(exist(a,b)) return;
      _data[a][b] = (byte *) calloc( size(a,b)/8 + 1, sizeof(byte));
    }

    // destruct the axis [a][b] 
    inline void destruct(vindex a, vindex b)
    {
      if (!exist(a,b)) return;
      free(_data[a][b]);
      _data[a][b] = nullptr;
    }


    // helper functions from parent
    using S::contains;
    using S::set;
    using S::reset;

    /**
     *  Necessary implementations
     */
    inline bool contains(vindex a, vindex b, vindex c) const { return exist(a,b) && _data[a][b][c/ 8] & (1 << (c % 8)); }
    inline void set(vindex a, vindex b, vindex c) { create(a,b);  _data[a][b][c/8] |= 1 << (c % 8); }
    inline void reset(vindex a, vindex b, vindex c) { create(a,b);  _data[a][b][c/8] &= ~(1 << (c%8)); }

    inline vindex size() const { return _size[0]; } 
    inline vindex size(vindex a) const { return _size[1] - a; }
    inline vindex size(vindex a, vindex b) const { return _size[2] - a - b; } 

    // helper functions for every set
    size_t count(vindex i, vindex j) const;
    bool empty(vindex i, vindex j) const;
    bool toFile(const std::string &filename, bool sparse = false) const;

    // delete all empty rows
    void compress();
};
