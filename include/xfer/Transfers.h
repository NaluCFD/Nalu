/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Transfers_h
#define Transfers_h

#include <Enums.h>

// yaml for parsing..
#include <yaml-cpp/yaml.h>
#include <xfer/Transfer.h>
#include <NaluParsing.h>

#include <map>
#include <string>
#include <vector>

namespace YAML {
class Node;
}

namespace sierra{
namespace nalu{

typedef std::vector<Transfer *> TransferVector;

class Simulation;

class Transfers : public TransferVector {
public:
  Transfers(Simulation& sim);
  ~Transfers();

  void load(const YAML::Node & node);
  void breadboard();
  void initialize();
  Simulation *root();
  Simulation *parent();

  Simulation &simulation_;

};

} // namespace nalu
} // namespace Sierra

#endif
