#include <boost/program_options.hpp>
#include <stdexcept>
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include "parallel.h"
#include "vector.h"
#include "domain.h"
#include "timer.h"
#include "tetrahedron.h"
#include "set.h"
#include "cubeset.h"
#include "squareset.h"
#include "filter.h"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  int scale;
  string tmpdir; //= "/local/rvveneti/";
  string findir; //= "/var/scratch/rvveneti/";
  string logdir; //= "log/";
  string file;
  string domain;
  bool sparse=false;


  // option description
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "display this message")
    ("scale,s", po::value<int>(&scale)->default_value(0), "set scale")
    ("domain,d", po::value<std::string>(&domain)->default_value("fund"), "domain")
    ("sparse,sp", "use sparse instead of full storage")
    ("file,f",  po::value<std::string>(&file)->default_value(""), "load from this file")
    ("legacy", "load using old file format")
    ("tdir",  po::value<std::string>(&tmpdir)->default_value("/local/rvveneti/"),
                "directory for temporary files")
    ("fdir",  po::value<std::string>(&findir)->default_value("/var/scratch/rvveneti/"),
               "directory for the final files")
    ("ldir",  po::value<std::string>(&logdir)->default_value("log/"),
                "directory for log files");


  po::positional_options_description p;
  p.add("scale", 1);

  // store the options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);


  // check if help is set
  if (vm.count("help")) {
    cout << desc<< endl;
    return 1;
  }


  if (scale ==0 && file.empty())
  {
    cout << desc << endl;
    cout << "Specify the scale for which we should gather or the filename from which we should continue" << endl;
    return 1;
  }

  // check if we should use sparse dataset instead
  sparse = (vm.count("sparse") > 0);

  // print settings
  cout << scale << endl << sparse << endl << file << endl << tmpdir << endl << findir << endl << logdir << endl << endl;

#ifdef USE_OMP
  #pragma omp parallel for schedule(static,1)
  for (int i = 0; i < 70; i++)
    printf("(%d,%d)", omp_get_thread_num(),omp_get_num_threads());
  cout << endl;
#endif

  mkdir(logdir.c_str(), 0777);
  mkdir(tmpdir.c_str(), 0777 );
  mkdir(findir.c_str(), 0777);


  FCubeTSetFull *FCubeFull = nullptr;
  FCubeTSetSparse *FCubeSparse = nullptr;

  // if load from file, do that first
  if (!file.empty()) 
  {
    cout << endl << endl << "Loading data file: " << file << endl <<endl;
    try {
      if (sparse)
        FCubeSparse = new FCubeTSetSparse(file,vm.count("legacy"));
      else
        FCubeFull = new FCubeTSetFull(file,vm.count("legacy"));
    } catch( std::runtime_error e) {
      cerr << "Error: " << e.what() << endl;
      return 0;
    }
    scale = FCubeFull ? FCubeFull->scale() : FCubeSparse->scale();
  }

  // Redirect output
  string filename = "";
  if (sparse)
    filename = to_string(scale) + ".FCubeTSetSparse";
  else
    filename = to_string(scale) + ".FCubeTSetFull";

  ofstream outfile(logdir + filename + ".log", ios::app);
  auto coutbuf = cout.rdbuf(outfile.rdbuf());

  if  (file.empty()) {
    cout << endl << endl << "Gathering results for scale = " << scale << endl<<endl << endl;
    double init_start = omp_get_wtime();
    if (vm.count("sparse"))
      FCubeSparse = new FCubeTSetSparse(scale,true);
    else
      FCubeFull = new FCubeTSetFull(scale,true);

    cout << "Initalisation took " << omp_get_wtime() - init_start << " seconds" << endl << endl;
  } else {
    cout << endl << endl << endl << "Loaded data from file: " << file << endl;
    cout << "Continue with data set for scale = " << scale << endl << endl;
  }

  if (FCubeFull) 
  {
    TriangleFilter<FCubeTSetFull> filter(*FCubeFull, findir + filename , tmpdir + filename + ".tmp", 60*60);
    filter.filter();
  }
  else
  {
    TriangleFilter<FCubeTSetSparse> filter(*FCubeSparse, findir + filename, tmpdir + filename + ".tmp", 60*60);
    filter.filter();
  }
  /*
  set->print();
  set->toFile("/tmp/bla1.fund",true);
  filter.sweep();
  set->print();
  set->toFile("/tmp/bla2.fund",true);
  */

  // remove objects
  if (FCubeFull)
    delete FCubeFull;
  else
    delete FCubeSparse;

  cout.rdbuf(coutbuf);
  return 0;
}
