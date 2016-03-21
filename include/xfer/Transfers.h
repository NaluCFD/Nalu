/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Transfers_h
#define Transfers_h

#include <Enums.h>

// yaml for parsing..
#include <yaml-cpp/yaml.h>

#include <map>
#include <string>
#include <vector>

namespace YAML {
class Node;
}

namespace sierra{
namespace nalu{

class Simulation;
class Transfer;

class Transfers {
public:
  Transfers(Simulation& sim);
  ~Transfers();

  void load(const YAML::Node & node);
  void breadboard();
  void initialize();
  void execute(); // general method to execute all xfers (as apposed to Realm)
  Simulation *root();
  Simulation *parent();

  Simulation &simulation_;
  std::vector<Transfer *> transferVector_;
};

} // namespace nalu
} // namespace Sierra

#endif
