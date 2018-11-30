/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MaterialProperty_h
#define MaterialProperty_h

#include "Enums.h"

// yaml for parsing..
#include <yaml-cpp/yaml.h>

// basic c++
#include <map>
#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class MaterialPropertys;
class MaterialPropertyData;
class ReferencePropertyData;
class PropertyEvaluator;

class MaterialProperty {
public:
  MaterialProperty(MaterialPropertys &materialPropertys, const std::string materialBlockName);
  
  ~MaterialProperty();
  
  void load(const YAML::Node & node);

  MaterialPropertys *parent();
  
  void extract_universal_constant( 
    const std::string name, double &value, const bool useDefault);
 
  MaterialPropertys &materialPropertys_;
  const std::string materialBlockName_;

  // start the parameters
  std::string propertyTableName_;
  
  // vectors and maps required to manage full set of options
  std::vector<std::string> targetNames_;
  std::map<std::string, double> universalConstantMap_;
  std::map<PropertyIdentifier, MaterialPropertyData*> propertyDataMap_;
  std::map<std::string, ReferencePropertyData*> referencePropertyDataMap_; /* defines overall species ordering */
  std::map<PropertyIdentifier, PropertyEvaluator*> propertyEvalMap_;
  std::map<std::string, ReferencePropertyData*> tablePropertyMap_;
};

} // namespace nalu
} // namespace Sierra

#endif
