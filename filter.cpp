#include <cstdio>
#include "parallel.h"
#include "vector.h"
#include "filter.h"
#include "domain.h"
#include "tetrahedron.h"
#include "cubeset.h"
#include "squareset.h"


// define conformity checker for triangle given by vertices (a,b,c)
template<typename T>
inline bool TriangleFilter<T>::valid(const Triangle<T::dim> &triangle) const
{
  // calculate edges of this triangle
  Simplex<3,3> edges = triangle.edges();
  if (!Triangle<3>::acute(edges)) return false;

  // calculate the normal
  Vector<3> normal = cross(edges[1],edges[0]);

  // calculate the constant specific for the triangle plane
  int d = dot(normal, triangle[0]);

  //Calculate the vector perpendicular to each side and the normal. Normals on each side of the prism
  Simplex<3,3> edge_normals;
  Vector<3> side_d;  //The constant expression for the plane equations of the sides
  //Convention, need to explain why it works?? Third point must be on the other side!!

  edge_normals = {{cross(normal, edges[0]),
                   cross(edges[1], normal),
                   cross(normal, edges[2])}};

  side_d ={{dot(edge_normals[0], triangle[0]),
            dot(edge_normals[1], triangle[0]),
            dot(edge_normals[2], triangle[2])}};

  /*
   * Initalize to zero
   */
  bool acute_above = false, acute_below = false;
  bool boundary = _domain.boundary(triangle);

  /*
   * In case we work store the acute indices, we need to locally store
   * if we found above and below. As we later put them into the data struct.
   */
  // count down, for higher performance
  for (int i = _domain.size() - 1; i >= 0; i--)
  {
    int dotprod = dot(_domain[i], normal);  //Calculate angle between normal and vector from plane to cube_pt
    //Check if current point lies on side of the triangle which isn't sharp yet and check if point lies in the
    //prism of possible sharp tetrahedron
    if (((dotprod > d && !acute_above) ||  //Dotprod > tri_d implies that the point lies "above" the triangle (same side as normal) 
          (dotprod < d && !acute_below)) && //Check if point lies in the prism
           (dot(_domain[i], edge_normals[0]) < side_d[0]) &&
           (dot(_domain[i], edge_normals[1]) < side_d[1]) &&
           (dot(_domain[i], edge_normals[2]) < side_d[2]) && //Dotprod < tri_d implies lies below
           (Tetrahedron::acute(triangle, _domain[i])) &&
           (_set.contains(triangle, _domain[i]))) // all the facets are in the set
    { //Passed all tests, we have found a correct tetrahedron on this side.
      if (dotprod > d) 
        acute_above = true;
      else 
        acute_below = true;
      if ((acute_above && acute_below) || boundary)
        return true;    
    }
  }
  return false;
}
template<>
inline bool TriangleFilter<SquareTSet>::valid(const Triangle<2> &triangle) const
{
  Simplex<2,3> edges = triangle.edges();
  if (!Triangle<2>::acute(edges)) return false;

   //Normals on the edges
  Vector<2> edge_normals[3]= {normal(edges[0]), normal(edges[1]), normal(edges[2])};

  // Line plane
  Vector<3> side_d ={{dot(edge_normals[0], triangle[0]),
                      dot(edge_normals[1], triangle[0]),
                      dot(edge_normals[2], triangle[2])}};
  
  Vector<3> orientation = {{dot(edge_normals[0], triangle[2] - triangle[0]),
                            dot(edge_normals[1], triangle[1] - triangle[0]),
                            dot(edge_normals[2], triangle[0] - triangle[2])}};

  // fix outgoing orientation
  for (size_t e = 0; e < 3; e++)
    if (orientation[e] > 0) {
      side_d[e] = -side_d[e];
      edge_normals[e] = Vector<2> {{0,0}} -edge_normals[e];
    }

  for (size_t e = 0; e < 3; e++) { // for each edge
    bool acute = false;
    Simplex<2,2> edge;
    if      (e == 0) edge = {{triangle[0], triangle[1]}};
    else if (e == 1) edge = {{triangle[0], triangle[2]}};
    else if (e == 2) edge = {{triangle[1], triangle[2]}};
    if  (_domain.boundary(edge)) continue;

    for (size_t i = 0; i < _domain.size(); i++) // for each point 
      if (dot(_domain[i], edge_normals[e]) > side_d[e] && // correct side 
          Triangle<2>{{edge[0], edge[1], _domain[i]}}.acute() && // acute
          _set.contains(edge[0], edge[1], _domain[i]))
      {
        acute = true;
        break;
      }
    if (!acute) return false;
  }
  return true;
}


// do one iteration of conform checking
template<typename T>
bool TriangleFilter<T>::sweep()
{
  bool changed = 0;
  double time_start = omp_get_wtime(), time_save = _interval;
  
#ifdef USE_OMP
  #pragma omp parallel for schedule(dynamic) 
#endif
  for (vindex i = 0; i < _set.size(); i++) {
    for (vindex j = 0; j < _set.size(i); j++) {
      for (vindex k = 0; k < _set.size(i,j); k++)
      {
        if (!_set(i,j,k)) continue;
        Triangle<3> triangle{{_set.vertex(i), _set.vertex(i+j),  _set.vertex(i+j+k)}};
        if (!valid(triangle)) { //remove from conf_mem_list
          changed = true;
          _set.reset(triangle);
        }
      }
      if (!_tmpfile.empty() && _interval > 0 &&  //Do we want to save the file?
            omp_get_thread_num() == 0 && //Only let master save to the file
          ((omp_get_wtime() - time_start) > time_save))  //Time to save current progress
      {  
        cout << "Saving: ";
        _set.print();
        if(!_set.toFile(_tmpfile + ".partial", true)) {
          cout << "Failed to save, try again at half the interval" << endl;
          time_save += _interval / 2;
        } else {
          rename( (_tmpfile + ".partial").c_str(), _tmpfile.c_str());
          time_save += _interval;
        }
      }
    }
  }
  return changed;
}

// do one iteration of cosy checking
template<>
bool TriangleFilter<SquareTSet>::sweep()
{
  bool changed = 0;
  const Square &square = _domain;
  const auto &set = _set;
#ifdef USE_OMP
  #pragma omp parallel for schedule(dynamic) 
#endif
  for (vindex i = 0; i < set.size(); i++) {
    for (vindex j = 0; j < set.size(i); j++) {
      for (vindex k = 0; k < set.size(i,j); k++)
      {
        if (!_set(i,j,k)) continue;
        Triangle<2> triangle{{square[i], square[i+j], square[i+j+k]}};
        if (!valid(triangle)) { //remove from conf_mem_list
          changed = true;
          _set.reset(i,j,k);
        }
      }
    }
  }
  return changed;
}

template<>
bool TriangleFilter<SquareTSet>::filter()
{
  cout << "\tFiltering " << _set.name() << " of scale " << (int) _domain._scale << " (p=" << _domain._scale - 1 << ")" << endl;
  double ftimer = omp_get_wtime();
  bool changed = true;
  size_t number = 0;
  while (changed) 
  {
    cout << "\t";
    number = _set.print();
    if (number == 0) break;

    double timer = omp_get_wtime();
    changed = sweep();
    cout << "\tSweep took " << omp_get_wtime() - timer << " seconds." <<  endl;
  }
  cout <<  "\tDone filtering. Took " << omp_get_wtime() - ftimer << " seconds." << endl;
  return (number > 0);
}

template<typename T>
inline bool TriangleFilter<T>::filter()
{
  cout << "Filtering " << _set.name() << " of scale " << (int) _domain._scale << " (p=" << _domain._scale - 1 << ")" << endl << endl;
  double ftimer = omp_get_wtime();
  size_t number = 0;
  bool changed = true;
  while (changed) 
  {
    number = _set.print();
    if (number == 0) break;
    cout << endl << "Gather and filter boundary facets for side 0" << endl;
    auto bdrtriangles = getboundaryfacets(0);
    TriangleFilter<decltype(bdrtriangles)> bdrfilter(bdrtriangles);
    bool removedbdr = bdrfilter.filter();

    if (removedbdr) {
      cout << endl << "\tRemove non-cosy boundary triangles" << endl; 
      setboundaryfacets(bdrtriangles, 0);
      cout << "\t";
      _set.print();
    }

    cout << endl << "Done filtering the boundaries, sweep entire set " << endl;
    double timer = omp_get_wtime();
    changed = sweep();
    cout << "Sweep took " << omp_get_wtime() - timer << " seconds." <<  endl;

    // compress the dataset, only works if this is implemented
    _set.compress();
    cout << "End of compressing" << endl;
  }
  cout << endl << endl <<  "Done filtering. Took " << omp_get_wtime() - ftimer << " seconds." << endl;
  _set.print();
  if (!_filename.empty()) _set.toFile(_filename, true);
  return (number > 0);
}

template<typename T>
void TriangleFilter<T>::setboundaryfacets(const SquareTSet &facets, byte axis, byte pt)
{
  if (facets.scale() != _set.scale()) return;

  const Square &square = facets.domain();

  // loop over every triangle in this square
  for(vindex i = 0; i < square.size(); i++)
    for (vindex j =i+1; j < square.size(); j++) 
      for(vindex k = j+1; k < square.size(); k++)
      {
        // if this triangle is not set in the facet set, remove the lifted version from the cube set
        if (!facets.contains(i,j-i,k-j)) 
          _set.reset(Triangle<3>{{square[i].lift(axis, pt), square[j].lift(axis,pt), square[k].lift(axis,pt)}});
      }
}

template<typename T>
SquareTSet TriangleFilter<T>::getboundaryfacets(byte axis, byte pt) const
{
  Square square(_set.scale());
  SquareTSet result(_set.scale(), false);

  // loop over every triangle in this square
  for(vindex i = 0; i < square.size(); i++)
    for (vindex j =i+1; j < square.size(); j++) 
      for(vindex k = j+1; k < square.size(); k++)
      {
        // we must now lift every point of the 2D triangle to the 3D triangle
        if (_set.contains(square[i].lift(axis, pt), square[j].lift(axis,pt), square[k].lift(axis,pt))) {
          //Triangle<2> {square[i], square[j],square[k]}.print();
          result.set(i, j-i, k-j);
        }
      }
  return result;
}

template class TriangleFilter<FCubeTSetFull>;
template class TriangleFilter<FCubeTSetSparse>;
