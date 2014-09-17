/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MaterialPropertys_h
#define MaterialPropertys_h

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
class MaterialProperty;
class MaterialPropertyData;
class ReferencePropertyData;
class PropertyEvaluator;
class Simulation;

typedef std::vector<MaterialProperty *> MaterialPropertyVector;

class MaterialPropertys : public MaterialPropertyVector {
public:
  MaterialPropertys(Realm& realm);
  
  ~MaterialPropertys();
  
  void load(const YAML::Node & node);
  
  void breadboard(){};
  
  Realm &realm_;
  std::string propertyTableName_;

  // vectors and maps required to manage full set of options
  std::vector<std::string> targetNames_;
  std::map<std::string, double> universalConstantMap_;
  std::map<PropertyIdentifier, MaterialPropertyData*> propertyDataMap_;
  std::map<std::string, ReferencePropertyData*> referencePropertyDataMap_; /* defines overall species ordering */
  std::map<PropertyIdentifier, PropertyEvaluator*> propertyEvalMap_;
  std::map<std::string, ReferencePropertyData*> tablePropertyMap_;
  Simulation *root();
  Realm *parent();  
  
};


} // namespace nalu
} // namespace Sierra

#endif
