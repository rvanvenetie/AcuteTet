#include <fstream>
#include <limits>
#include "set.h"

#define tempTD template <typename T, typename D>

// initalize
tempTD void TriangleSet<T,D>::init(bool set)
{
  _data = (byte ***) calloc(size(), sizeof(byte **));
  //_data2 = (boost::dynamic_bitset<> **) calloc(size(), sizeof(boost::dynamic_bitset<> *));
  for (vindex i = 0; i < size(); i++) {
    _data[i] = (byte **) calloc( size(i), sizeof (byte *));
    //_data2[i] = (boost::dynamic_bitset<> *) calloc( size(i), sizeof (boost::dynamic_bitset<>));
    for (vindex j =0; j < size(i); j++) {
      // calloc, bits
      _data[i][j] = (byte *) calloc( size(i,j)/8 + 1, sizeof(byte));
      //new(&_data2[i][j]) boost::dynamic_bitset<>(size(i,j));
      if (set)
        for (vindex k =0; k < size(i,j); k++) 
          this->set(i,j,k);
    }
  }
}
// destructor
tempTD TriangleSet<T,D>::~TriangleSet()
{
  for (vindex i = 0;i < size(); i++) {
    for (vindex j =0; j < size(i); j++)  {
      //_data2[i][j].~dynamic_bitset<>();
      free(_data[i][j]);
    }
    free(_data[i]);
    //free(_data2[i]);
  }
  free(_data);
  //free(_data2);
}

// count number of triangles
tempTD size_t TriangleSet<T,D>::count(vindex i, vindex j) const
{
  size_t result = 0;
  for (vindex k =0; k < size(i,j) / 8 + 1; k++)
    result += std::bitset<8>(_data[i][j][k]).count();
  return result;
}

tempTD size_t TriangleSet<T,D>::count() const
{
  size_t result = 0;
  for (vindex i = 0; i < size(); i++) 
    for (vindex j =0; j < size(i); j++) 
      result += count(i,j);
  return result;
}

// count memory
tempTD size_t TriangleSet<T,D>::memory() const
{
  size_t result = 0;
  result += size() * sizeof(byte **); // first dimension
  for (vindex i = 0; i <size(); i++) 
  {
    result += size(i) * sizeof(byte *);
    for (vindex j =0 ; j < size(i); j++)
      result += (size(i,j)/ 8 + 1) * sizeof(byte);
  }
  return result;  
}

// debug
tempTD size_t TriangleSet<T,D>::print() const
{ 
  std::string units[] = {"B", "kB", "MB", "GB"};
  size_t count = 0;
  vindex size = memory();
  vindex i = 0;
  while (size > 1024) {
    size /= 1024;
    i++;
  }
  cout << _name << " contains " <<  (count = this->count()) << " triangles. Memory used is " << size << units[i] << endl;
  return count;
}

// to file
tempTD bool TriangleSet<T,D>::toFile(const std::string &filename, TriangleSet<T,D>::FileOptions options) const 
{
  ofstream file(filename, ios::trunc  | ios::binary);

  if (!file.good()) return false;

  // Write the name and scale to the file
  file << _name << " " << (int) _scale << '\n';

  //Saves to file by looping over first two axis.
  //Then it saves <c><last_dimension>
  //Here c indicates whether the row is saved (0 can occur if MEM_LIST_SAVE_CLEAN and row is empty)
  for (vindex i = 0;i < size(); i++) {
    for (vindex j =0; j < size(i); j++)  {
      if (options == FULL || (options == SPARSE && count(i,j))) 
      {
        file << (byte) 1;
        file.write((const char *) _data[i][j],size(i,j)/8 + 1);
        if (!file.good()) return false;
      } else {
        file << (byte) 0;
      }
    }
  }
  return true;
}

// from file
tempTD bool TriangleSet<T,D>::fromFile(ifstream &file)
{
  if (!file.good()) return false;
  for (vindex i = 0; i < size(); i++) {
    for (vindex j =0; j < size(i); j++) {
      byte last;
      file >> last;
      if (last == 1)
      {
        file.read((char *) _data[i][j], size(i,j)/8 + 1);
        if (!file.good()) return false;
      }
    }
  }
  return  true;
}

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
  for (vindex cidx =0 ; cidx < _domain.size(); cidx++) {
    // so cidx maps to fidx
    vindex fidx = _index[cidx];

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
    for (vindex idx =0; idx < _domain.size(); idx++) {
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
/**
 *  Explicit initalisation of the template classes
 */
template class TriangleSet<CubeTriangleSet, Cube>;
template class TriangleSet<FundcubeTriangleSet, Cube>;
