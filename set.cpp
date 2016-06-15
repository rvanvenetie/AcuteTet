#include <fstream>
#include <limits>
#include "set.h"
#include "cubeset.h"
#include "squareset.h"

#define tempTD template <typename T, typename D>

// initalize
tempTD void TriangleSet<T,D>::init(bool set)
{
  _data = (byte ***) calloc(size(), sizeof(byte **));
  //_data2 = (boost::dynamic_bitset<> **) calloc(size(), sizeof(boost::dynamic_bitset<> *));
  for (size_t i = 0; i < size(); i++) {
    _data[i] = (byte **) calloc( size(i), sizeof (byte *));
    //_data2[i] = (boost::dynamic_bitset<> *) calloc( size(i), sizeof (boost::dynamic_bitset<>));
    for (size_t j =0; j < size(i); j++) {
      // calloc, bits
      _data[i][j] = (byte *) calloc( size(i,j)/8 + 1, sizeof(byte));
      //new(&_data2[i][j]) boost::dynamic_bitset<>(size(i,j));
      if (set)
        for (size_t k =0; k < size(i,j); k++) 
          this->set(i,j,k);
    }
  }
}
// destructor
tempTD TriangleSet<T,D>::~TriangleSet()
{
  for (size_t i = 0;i < size(); i++) {
    for (size_t j =0; j < size(i); j++)  {
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
  for (size_t i = 0; i < size(); i++) 
    for (size_t j =0; j < size(i); j++) 
      result += count(i,j);
  return result;
}

// count memory
tempTD size_t TriangleSet<T,D>::memory() const
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
tempTD size_t TriangleSet<T,D>::print() const
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
  for (size_t i = 0;i < size(); i++) {
    for (size_t j =0; j < size(i); j++)  {
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
  for (size_t i = 0; i < size(); i++) {
    for (size_t j =0; j < size(i); j++) {
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

/**
 *  Explicit initalisation of the template classes
 */
template class TriangleSet<CubeTriangleSet, Cube>;
template class TriangleSet<FundcubeTriangleSet, Cube>;
template class TriangleSet<SquareTriangleSet, Square>;
