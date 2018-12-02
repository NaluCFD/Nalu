/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MaterialPropertys_h
#define MaterialPropertys_h

#include "Enums.h"

// yaml for parsing..
#include <yaml-cpp/yaml.h>

// basic c++
#include <vector>

namespace YAML {
class Node;
}

namespace sierra{
namespace nalu{

class Realm;
class MaterialProperty;
class Simulation;

typedef std::vector<MaterialProperty *> MaterialPropertyVector;

class MaterialPropertys {
public:
  MaterialPropertys(Realm& realm);
  
  ~MaterialPropertys();
  
  void load(const YAML::Node & node);
  
  void breadboard(){};
  
  // ease of access methods to particular initial condition
  size_t size() {return materialPropertyVector_.size();}
  MaterialProperty *operator[](int i) { return materialPropertyVector_[i];}

  Simulation *root();
  Realm *parent();  

  Realm &realm_;

  MaterialPropertyVector materialPropertyVector_;
  std::vector<std::string> targetNames_;
};


} // namespace nalu
} // namespace Sierra

#endif
