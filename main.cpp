#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <omp.h>
#include "vector.h"
#include "domain.h"
#include "timer.h"
#include "tetrahedron.h"
#include "set.h"
#include "cubeset.h"
#include "squareset.h"
#include "filter.h"
using namespace std;

int main(int argc, char *argv[]) {
  std::string tmpdir = "/local/rvveneti/";
  std::string findir = "/var/scratch/rvveneti/";
  std::string logdir = "log/";
  cout << argc << endl;
  if (!(argc == 4 || argc == 2))
  {
    cout << "Illegal parameters. " << endl;
    cout << "\t <scale>" << endl;
    cout << "\t <filename>" << endl;
    cout << "\t <scale> <tmpdir> <finaldir>" << endl;
    cout << "\t <filename> <tmpdir> <finaldir>" << endl;
    return 0;
  }
  if (argc == 4) {
    tmpdir = argv[2];
    findir = argv[3];
  } 
  string loadfilename;
  int scale = 0;
  if (isdigit(argv[1][0]))
    scale = atoi(argv[1]);
  else
    loadfilename = argv[1];

  string filename = "fund_" + string(argv[1]);


  #pragma omp parallel for schedule(static,1)
  for (int i = 0; i < 70; i++)
    printf("(%d,%d)", omp_get_thread_num(),omp_get_num_threads());
  cout << endl;

  mkdir(logdir.c_str(), 0777);
  mkdir(tmpdir.c_str(), 0777 );
  mkdir(findir.c_str() , 0777);


  // Redirect output
  //ofstream outfile(logdir + filename + ".log", ios::app);
  //auto coutbuf = cout.rdbuf(outfile.rdbuf());
  FCubeTSet<TFullSet<3>> *set;
  if  (scale) {
    cout << endl << endl << "Gathering results for scale = " << scale << endl<<endl << endl;
    double init_start = omp_get_wtime();
    set = new FCubeTSet<TFullSet<3>>(scale,true);
    cout << "Initalisation took " << omp_get_wtime() - init_start << " seconds" << endl << endl;
  } else {
    cout << endl << endl << "Loading data file: " << loadfilename << endl <<endl;
    set = new FCubeTSet<TFullSet<3>>(loadfilename, true);
  }

  TriangleFilter<FCubeTSet<TFullSet<3>>> filter(*set, findir + filename +".fund" , tmpdir + filename + ".tmp.fund", 60*60);
  filter.filter();

  // remove objects
  delete set;
  //cout.rdbuf(coutbuf);
  return 0;
}
