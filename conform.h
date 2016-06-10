#pragma once
#include <string>

/**
 *  Base class for conform set generator
 */
class Conform {
  private:
    std::string _filename;
  public:
};

/**
 *  Class that generates the conform set for the unit cube.
 */
class ConformCube : public Conform {

  public:
    void fund();
};
