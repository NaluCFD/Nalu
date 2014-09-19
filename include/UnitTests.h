/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef UnitTests_h
#define UnitTests_h

#include <Simulation.h>

namespace sierra{
namespace nalu{

class UnitTests {
public:
  UnitTests(Simulation& sim) : sim_(sim) {}

  ~UnitTests(){}

  void load(const YAML::Node & node);
  //void breadboard();
  //void initialize();
  void run();
  Simulation *root() { return &sim_; }
  Simulation *parent() { return &sim_; }

  Simulation& sim_;
};

} // namespace nalu
} // namespace Sierra

#endif

