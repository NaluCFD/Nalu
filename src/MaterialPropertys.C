/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "MaterialPropertys.h"

#include "MaterialProperty.h"
#include "Realm.h"

// yaml for parsing..
#include <yaml-cpp/yaml.h>

// basic c++
#include <vector>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// MaterialPropertys - manager of all material property types
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MaterialPropertys::MaterialPropertys(Realm& realm)
  : realm_(realm)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
MaterialPropertys::~MaterialPropertys()
{ 
  for ( size_t i = 0; i < materialPropertyVector_.size(); ++i ) {
    delete materialPropertyVector_[i];
  }
}

//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
MaterialPropertys::load(const YAML::Node & node) 
{
  const YAML::Node y_material_propertys = node["material_properties"];
  if (y_material_propertys) {
    
    // support two input file options
    const YAML::Node y_props = expect_sequence(y_material_propertys, "properties", true);
    if (y_props) {
      for (size_t iprop = 0; iprop < y_props.size(); ++iprop) {
        const YAML::Node y_prop = y_props[iprop] ;

        // extract name
        std::string materialBlockName = "na";
        get_required(y_prop, "name", materialBlockName);
          
        // create a new material property object
        MaterialProperty *matPropBlock = new MaterialProperty(*this, materialBlockName);
        
        // load and push back
        matPropBlock->load(y_prop);
        materialPropertyVector_.push_back(matPropBlock);
        
        // push back target names to material propertys
        targetNames_.insert(targetNames_.begin(), matPropBlock->targetNames_.begin(), matPropBlock->targetNames_.end()); 
      }
    }
    else {
      // create a new material property object with default name
      MaterialProperty *matPropBlock = new MaterialProperty(*this, "default");
      
      // load and push back
      matPropBlock->load(y_material_propertys);
      materialPropertyVector_.push_back(matPropBlock);
      
      // push back target names to material propertys
      targetNames_.insert(targetNames_.begin(), matPropBlock->targetNames_.begin(), matPropBlock->targetNames_.end()); 
    }
  }
  else {
    throw std::runtime_error("Error: material_properties::load(): 'material_properties' line does not exist");
  }
}

Simulation* MaterialPropertys::root() { return parent()->root(); }
Realm *MaterialPropertys::parent() { return &realm_; }

} // namespace nalu
} // namespace Sierra
