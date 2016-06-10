#include <iostream>
#include "vector.h"
#include "timer.h"
#include "tetrahedron.h"
using namespace std;

int main() {
  Vector<3> a{1,2,3};
  Vector<3> b{0,-5,0};
  Timer tm;
  int j;
  srand(25);
  j = 0;
  tm.start();
  for (int i = 0; i < 500000; i++)
  {
    Tetrahedron c = Tetrahedron::random(25);
    if (c.acute()) j++;
  }
  tm.stop();
  std::cout << j << " " << tm.duration() << endl;
  a.print();
  std::cout << "LOLOLOL" << std::endl;
  std::cout << sizeof(a) << std::endl;
  std::cout << a.zero() << std::endl;
  std::cout << b.zero() << std::endl;
  std::cout << (a == b) << std::endl;
  std::cout << dot(a,b) << std::endl;
  std::cout << a.dot(b) << std::endl;
  (a*b).print();
  //c.print();
  //std::cout << c.acute() << std::endl;

  return 0;
}
