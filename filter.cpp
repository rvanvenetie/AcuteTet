#include <omp.h>
#include "vector.h"
#include "filter.h"
#include "domain.h"
#include "tetrahedron.h"


// define conformity checker for triangle given by vertices (a,b,c)
template<typename T>
inline bool TriangleFilter<T>::conform(const Triangle<3> &triangle) const
{
  // calculate edges of this triangle
  Polygon<3,3> edges = triangle.edges();
  if (!triangle.acute(edges)) return false;

  // calculate the normal
  Vector<3> normal = cross(edges[1],edges[0]);

  // calculate the constant specific for the triangle plane
  int d = dot(normal, triangle[0]);

  //Calculate the vector perpendicular to each side and the normal. Normals on each side of the prism
  Polygon<3,3> edge_normals;
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
  for (size_t i = 0; i < _domain.size(); i++)
  {
    int dotprod = dot(_domain[i], normal);  //Calculate angle between normal and vector from plane to cube_pt
    //Check if current point lies on side of the triangle which isn't sharp yet and check if point lies in the
    //prism of possible sharp tetrahedron
    if (((dotprod > d && !acute_above) ||  //Dotprod > tri_d implies that the point lies "above" the triangle (same side as normal) 
          (dotprod < d && !acute_below)) && //Check if point lies in the prism
           dot(_domain[i], edge_normals[0]) < side_d[0] &&
           dot(_domain[i], edge_normals[1]) < side_d[1] &&
           dot(_domain[i], edge_normals[2]) < side_d[2] && //Dotprod < tri_d implies lies below
           Tetrahedron::acute(triangle, _domain[i]) &&
          _set.contains(triangle, _domain[i])) // all the facets are in the set
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

// do one iteration of conform checking
template<>
bool TriangleFilter<CubeTriangleSet>::sweep() 
{
  bool changed = 0;
  const Cube &cube = _domain;
  for (size_t i = 0; i < cube.size(); i++) 
    for (size_t j = 0; j < cube.size() - i; j++)
      for (size_t k = 0; k < cube.size() - j - i; k++)
      {
        if (!_set(i,j,k)) continue;
        if (!conform(Triangle<3>{cube[i], cube[i+j],  cube[i+j+k]})) { //remove from conf_mem_list
          changed = true;
          _set.reset(i,j,k);
        }

      }
  return changed;
}

// do one iteration of conform checking
template<>
bool TriangleFilter<FundcubeTriangleSet>::sweep()
{
  bool changed = 0;
  const Cube &cube = _domain;
  const Fundcube &fund = _set.fund();
  double time_start = omp_get_wtime(), time_save = _interval;
  

  #pragma omp parallel for schedule(dynamic) 
  for (vindex i = 0; i < fund.size(); i++) {
    for (vindex j = 0; j < cube.size() - i; j++) {
      for (vindex k = 0; k < cube.size() - j - i; k++)
      {
        if (!_set(i,j,k)) continue;
        Triangle<3> triangle{_set.vertex(i), _set.vertex(i+j),  _set.vertex(i+j+k)};
        if (!conform(triangle)) { //remove from conf_mem_list
          changed = true;
          _set.reset(triangle);
        }
      }
      if (!_tmpfile.empty() && _interval > 0 &&  //Do we want to save the file?
            omp_get_thread_num() == 0 && //Only let master save to the file
          ((omp_get_wtime() - time_start) > time_save))  //Time to save current progress
      {  
        cout << "Saving tmp file. ";
        _set.print();
        if(!_set.toFile(_tmpfile, TriangleSet<FundcubeTriangleSet, Cube>::SPARSE)) {
          cout << "Failed to save, try again at half the interval" << endl;
          time_save += _interval / 2;
        } else 
          time_save += _interval;
      }
    }
  }
  return changed;
}

template<typename T>
inline bool TriangleFilter<T>::filter()
{
  cout << "Filtering a triangle set of scale " << (int) _domain._scale << endl;
  double ftimer = omp_get_wtime();
  size_t number = 0;
  bool changed = true;
  while (changed) 
  {
    number = _set.print();
    cout << endl;
    double timer = omp_get_wtime();
    changed = sweep();
    cout << endl << "Sweep took " << omp_get_wtime() - timer << " seconds." <<  endl;
  }
  _set.print();
  cout << endl <<  "Done filtering. Took " << omp_get_wtime() - ftimer << " seconds." << endl << endl;
  _set.toFile(_filename);
  return (number > 0);
}

template class TriangleFilter<CubeTriangleSet>;
template class TriangleFilter<FundcubeTriangleSet>;
