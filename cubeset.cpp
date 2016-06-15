#include <fstream>
#include <limits>
#include "cubeset.h"

CubeTriangleSet * CubeTriangleSet::fromFile(const std::string &filename) 
{
  std::string classname;
  int scale;

  ifstream file(filename, ios::binary);
  if (!file.good()) return nullptr;
  file >> classname >>  scale;
  // read until end of line
  file.ignore(123456, '\n');
  if (classname != "CubeTriangleSet") return nullptr;
  CubeTriangleSet * result = new CubeTriangleSet(scale,false);
  if (result->TriangleSet<CubeTriangleSet, Cube>::fromFile(file))
    return result;
  else
    return nullptr;
}

FundcubeTriangleSet::FundcubeTriangleSet(byte scale, bool set) :
      TriangleSet<FundcubeTriangleSet, Cube>(scale),
      _fund(scale),
      _vertices(_domain.size()),
      _index(_domain.size(),numeric_limits<vindex>::max()),
      _symindex(_domain.size()),
      _sympt(_domain.size()*48)
{ 
  _cubesize = _domain.size();
  _fundsize = _fund.size();
  this->_name = "FundcubeTriangleSet";
  this->init(set); 

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

FundcubeTriangleSet * FundcubeTriangleSet::fromFile(const std::string &filename) 
{
  std::string classname;
  int scale;

  ifstream file(filename, ios::binary);
  if (!file.good()) return nullptr;
  file >> classname >> scale;
  // read until end of line
  file.ignore(123456, '\n');
  if (classname != "FundcubeTriangleSet") return nullptr;
  FundcubeTriangleSet * result = new FundcubeTriangleSet(scale,false);
  if (result->TriangleSet<FundcubeTriangleSet, Cube>::fromFile(file))
    return result;
  else
    return nullptr;
}
