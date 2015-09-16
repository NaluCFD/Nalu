/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef InitialConditions_h
#define InitialConditions_h

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
class InitialConditions;
class Simulation;

class InitialCondition {
 public:
 InitialCondition(InitialConditions& ics) : initialConditions_(ics), theIcType_(UserDataType_END) {}
  
  virtual ~InitialCondition() {}
  
  InitialCondition * load(const YAML::Node & node) ;
  Simulation *root();
  InitialConditions *parent();
  
  void breadboard()
  {
    // nothing
  }
  
  InitialConditions& initialConditions_;
  
  std::string icName_;
  std::vector<std::string> targetNames_;
  UserDataType theIcType_;
};
 
 typedef std::vector<InitialCondition *> InitialConditionVector;
 
 class InitialConditions {
 public:
 InitialConditions(Realm& realm) : realm_(realm) {}
 
 ~InitialConditions() 
   {
     for ( size_t j_initial_condition = 0; j_initial_condition < initialConditionVector_.size(); ++j_initial_condition ) {
       delete initialConditionVector_[j_initial_condition];
     }
   }
 
 InitialConditions* load(const YAML::Node & node) 
 {
   InitialCondition tmp_initial_condition(*this);
   
   const YAML::Node *initial_conditions = node.FindValue("initial_conditions");
   if (initial_conditions) {
     for ( size_t j_initial_condition = 0; j_initial_condition < initial_conditions->size(); ++j_initial_condition ) {
       const YAML::Node & initial_condition_node = (*initial_conditions)[j_initial_condition];
       InitialCondition* ic = tmp_initial_condition.load(initial_condition_node);
       initialConditionVector_.push_back(ic);
     }
   }
   else
     throw std::runtime_error("parser error InitialConditions::load");
   
   return this;
 }
 
 void breadboard()
 {
   for ( size_t j_initial_condition = 0; j_initial_condition < initialConditionVector_.size(); ++j_initial_condition ) {
     initialConditionVector_[j_initial_condition]->breadboard();
   }
 }
 
 Simulation *root();
 Realm *parent();
 
 // ease of access methods to particular initial condition
 size_t size() {return initialConditionVector_.size();}
 InitialCondition *operator[](int i) { return initialConditionVector_[i];}
 
 Realm &realm_;
 InitialConditionVector initialConditionVector_;
}; 
 
} // namespace nalu
} // namespace Sierra

#endif
