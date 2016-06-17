#include <fstream>
#include <limits>
#include "set.h"
#include "cubeset.h"
#include "squareset.h"

#define tempTD template <typename T, byte D>
#define tempD  template <byte D>


// count number of triangles
tempTD size_t TSet<T,D>::count() const
{
  size_t result = 0;
  for (size_t i = 0; i < size(); i++) 
    for (size_t j =0; j < size(i); j++) 
      result += count(i,j);
  return result;
}

// count memory
tempTD size_t TSet<T,D>::memory() const
{
  size_t result = 0;
  result += size() * sizeof(byte **); // first dimension
  for (size_t i = 0; i <size(); i++) 
  {
    result += size(i) * sizeof(byte *);
    for (size_t j =0 ; j < size(i); j++)
      result += (size(i,j)/ 8 + 1) * sizeof(byte);
  }
  return result;  
}

// debug
tempTD size_t TSet<T,D>::print() const
{ 
  std::string units[] = {"B", "kB", "MB", "GB"};
  size_t count = 0;
  size_t size = memory();
  size_t i = 0;
  while (size > 1024) {
    size /= 1024;
    i++;
  }
  cout << _name << " contains " <<  (count = this->count()) << " triangles. Memory used is " << size << units[i] << endl;
  return count;
}

/**
 *  TFullSet implementation
 */
// initalization method
tempD void TFullSet<D>::init(tindex axis, byte scale, bool set) {
  _name = "TFullSet";
  _scale = scale;
  _size = axis;

  _data = (byte ***) calloc(size(), sizeof(byte **));
  for (size_t i = 0; i < size(); i++) {
    _data[i] = (byte **) calloc( size(i), sizeof (byte *));
    for (size_t j =0; j < size(i); j++) {
      // calloc, bits
      _data[i][j] = (byte *) calloc( size(i,j)/8 + 1, sizeof(byte));
      if (set)
        for (size_t k =0; k < size(i,j); k++) 
          this->set(i,j,k);
    }
  }
}

// from file
tempD void TFullSet<D>::init(ifstream &file, tindex axis, byte scale) 
{
  init(axis, scale, false);
  if (!file.good()) throw std::runtime_error("Failed to load from file");
  for (size_t i = 0; i < size(); i++) {
    for (size_t j =0; j < size(i); j++) {
      byte last;
      file >> last;
      if (last == 1)
      {
        file.read((char *) _data[i][j], size(i,j)/8 + 1);
        if (!file.good()) throw std::runtime_error("Failed to load from file");
      }
    }
  }
}

// destructor
tempD TFullSet<D>::~TFullSet()
{
  for (size_t i = 0;i <size(); i++) {
    for (size_t j =0; j < size(i); j++)  {
      free(_data[i][j]);
    }
    free(_data[i]);
  }
  free(_data);
}

// count the number of triangles for index i,j
tempD size_t TFullSet<D>::count(vindex i, vindex j) const
{
  size_t result = 0;
  for (vindex k =0; k < (_size[2]-i-j) / 8 + 1; k++)
    result += std::bitset<8>(_data[i][j][k]).count();
  return result;
}

// returns whether the axis is empty
tempD bool TFullSet<D>::empty(vindex i, vindex j) const
{
  for (vindex k =0; k < (_size[2]-i-j) / 8 + 1; k++)
    if(std::bitset<8>(_data[i][j][k]).count()) return false;
  return true;
}

// to file
tempD bool TFullSet<D>::toFile(const std::string &filename, bool sparse) const 
{
  ofstream file(filename, ios::trunc  | ios::binary);

  if (!file.good()) return false;

  // Write the name and scale to the file
  file <<  _name << " " << (int) _scale << " " << D << '\n';

  //Saves to file by looping over first two axis.
  //Then it saves <c><last_dimension>
  //Here c indicates whether the row is saved (0 can occur if MEM_LIST_SAVE_CLEAN and row is empty)
  for (size_t i = 0;i < size(); i++) {
    for (size_t j =0; j < size(i); j++)  {
      if (!sparse || (sparse && empty(i,j))) 
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

/**
 *  Explicit initalisation of the template classes
 */
template class TFullSet<3>;
template class TSet<TFullSet<3>,3>;
template class TFullSet<2>;
template class TSet<TFullSet<2>,2>;
