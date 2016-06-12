#include <iostream>
#include <vector>
#include <omp.h>
#include "vector.h"
#include "domain.h"
#include "timer.h"
#include "tetrahedron.h"
#include "set.h"
#include "filter.h"
using namespace std;

int main() {
  int scale = 7;
  cin >> scale;
  cout << endl << "Gathering results for scale = " << scale << endl<<endl;
  double init_start = omp_get_wtime();
  FundcubeTriangleSet set(scale,true);
  cout << "Initalisation took " << omp_get_wtime() - init_start << " seconds" << endl << endl;

  TriangleFilter<FundcubeTriangleSet> filter(set, "/tmp/final.fund", "/tmp/tst.fund", 50*50);
  filter.filter();
  return 0;
}
