/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef BoundaryConditions_h
#define BoundaryConditions_h

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

class Realm;
class BoundaryConditions;
class Simulation;

class BoundaryCondition {
 public:
 BoundaryCondition(BoundaryConditions& bcs) : boundaryConditions_(bcs) {}
  
  virtual ~BoundaryCondition() {}
  
  BoundaryCondition * load(const YAML::Node & node) ;
  Simulation *root();
  BoundaryConditions *parent();
  
  void breadboard()
  {
    // nothing
  }
  
  std::string bcName_;
  std::string targetName_;
  BoundaryConditionType theBcType_;
  BoundaryConditions& boundaryConditions_;
};
 
 typedef std::vector<BoundaryCondition *> BoundaryConditionVector;
 
 class BoundaryConditions {
 public:
   
 BoundaryConditions(Realm& realm) 
   : realm_(realm) {}
 ~BoundaryConditions() {
   for ( size_t iboundary_condition = 0; iboundary_condition < boundaryConditionVector_.size(); ++iboundary_condition ) {
     delete boundaryConditionVector_[iboundary_condition];
   }
 }

 BoundaryConditions* load(const YAML::Node & node) 
 {
   BoundaryCondition tmp_boundary_condition(*this);
   
   const YAML::Node *boundary_conditions = node.FindValue("boundary_conditions");
   if (boundary_conditions) {
     for ( size_t iboundary_condition = 0; iboundary_condition < boundary_conditions->size(); ++iboundary_condition ) {
       const YAML::Node & boundary_condition_node = (*boundary_conditions)[iboundary_condition];
       BoundaryCondition* bc = tmp_boundary_condition.load(boundary_condition_node);
       boundaryConditionVector_.push_back(bc);
     }
   }
   else
     throw std::runtime_error("parser error BoundaryConditions::load");
   
   return this;
 }
 
 void breadboard()
 {
   for ( size_t iboundary_condition = 0; iboundary_condition < boundaryConditionVector_.size(); ++iboundary_condition ) {
     boundaryConditionVector_[iboundary_condition]->breadboard();
   }
 }
 
 Simulation *root();
 Realm *parent();
 
 // ease of access methods to particular boundary condition
 size_t size() {return boundaryConditionVector_.size();}
 BoundaryCondition *operator[](int i) { return boundaryConditionVector_[i];}
 
 Realm &realm_;
 BoundaryConditionVector boundaryConditionVector_;
};

} // namespace nalu
} // namespace Sierra

#endif
