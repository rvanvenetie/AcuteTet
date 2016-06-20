#include <stdexcept>
#include <fstream>
#include <limits>
#include "cubeset.h"

// S is the implementor of the set
#define tempS template <typename S>
/**
 *  Definite CubeTSet
 */
// Constructor from file
tempS CubeTSet<S>::CubeTSet(std::string filename) : _domain(0)
{
  std::string classname;
  int scale;

  ifstream file(filename, ios::binary);
  if (!file.good()) throw std::runtime_error("Failed to open file");
  file >> classname >>  scale;
  // read until end of line
  file.ignore(123456, '\n');
  if (classname != "CubeTSet") throw std::runtime_error("File contains invalid class");
  _domain = Cube(scale);
  S::init({_domain.size(), _domain.size(), _domain.size()});
  S::_name = "CubeTSet";
}

/**
 *  Define FCubeTSet
 */
// this initalizes all the variables
tempS void FCubeTSet<S>::init(byte scale) 
{
  _cubesize = _domain.size();
  _fundsize = _fund.size();
  _vertices.reserve(_cubesize);
  _index.assign(_cubesize, numeric_limits<vindex>::max());
  _symindex.reserve(_cubesize);
  _sympt.reserve(_cubesize * 48);

  int cnt = 0;
  // create index array, first the fundamental domain
  for (auto const &vtx : _fund._vertices)
    _index[_domain.index(vtx)] = cnt++;

  for (auto const &vtx : _domain._vertices)
    if (_index[_domain.index(vtx)] == numeric_limits<vindex>::max())
      _index[_domain.index(vtx)] = cnt++;
  
  // do the reverse, translate fund idx -> cube points
  for (size_t cidx =0 ; cidx < _domain.size(); cidx++) {
    // so cidx maps to fidx
    size_t fidx = _index[cidx];

    _vertices[fidx] = _domain[cidx];
  }

  /*
   * For each vertex we give the symmetry number needed to transform this
   * point into the fundamental domain
   */
  for (auto const &vtx : _domain._vertices)
    _symindex[_domain.index(vtx)] = _fund.symindex(vtx);

  // for each symmetry and each point we generate the symmetry index
  for (size_t sym = 0; sym < 48; sym++)
    for (size_t idx =0; idx < _domain.size(); idx++) {
      _sympt[idx +sym*_domain.size()] = _index[_domain.index(_fund.symmetry(_domain[idx], sym))];
    }
}

tempS FCubeTSet<S>::FCubeTSet(byte scale, bool set) : _domain(scale), _fund(scale)
{ 
  init(scale);
  S::init({{(vindex)_fundsize, (vindex) _cubesize, _cubesize}}, scale, set);
  S::_name = "FCubeTSet";
}

tempS FCubeTSet<S>::FCubeTSet(const std::string &filename, bool legacy) :
  _domain(0), _fund(0)
{
  if (legacy) { legacyload(filename); return; }

  std::string classname;
  int scale;

  ifstream file(filename, ios::binary);
  if (!file.good()) throw std::runtime_error("Failed opening file");
  file >> classname >> scale;
  // read until end of line
  file.ignore(123456, '\n');
  if (classname != "FCubeTSet" && classname != "FundcubeTriangleSet")
    throw std::runtime_error("File contains wrong classtype");
  _domain = Cube(scale);
  _fund   = Fundcube(scale);
  // init the values
  init(scale); 
  // load from file
  S::init(file, {_fundsize, _cubesize, _cubesize}, scale); 
  S::_name = "FCubeTSet";
}


tempS void FCubeTSet<S>::legacyload(const std::string &filename)
{
  ifstream file(filename, ios::binary);
  int scale;
  if (!file.good()) throw std::runtime_error("Failed to load");
  char buffer[64];
  file.read (buffer,64);
  scale = buffer[56]+1;

  _domain = Cube(scale);
  _fund   = Fundcube(scale);
  // init the values
  init(scale); 
  // create the axis
  S::init( {_fundsize, _cubesize, _cubesize},scale, false); 

  byte * tmp = (byte *) calloc(S::size(0), sizeof(byte));
  size_t cnt = 0;
  for (size_t i = 0; i < S::size(); i++) {
    for (int j = 0; j < S::size(i)-1; j++) {
      byte last;
      file >> last;
      if (last == 1)
      {
        file.read((char *) tmp, (S::size(i,j)-2)/8 + 1);
        if (!file.good()) throw std::runtime_error("Failed to load legacy");
        for (int k=0; k < S::size(i,j)-2; k++) {
          if (tmp[k/8] & (1 << (k % 8))) {
            cnt++;
            S::set(i,j+1,k+1);
          }
        }
      }
    }
  }
  cout << "Loaded legacy for scale " << scale << " (p= " << scale - 1 << "). Amount of tri: " << cnt << endl << endl;
  free(tmp);
  S::_name = "FCubeTSet";
}

template class FCubeTSet<TFullSet<3>>;
