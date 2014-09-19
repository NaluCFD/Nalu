/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MaterialProperty_h
#define MaterialProperty_h

#include <Enums.h>

// yaml for parsing..
#include <yaml-cpp/yaml.h>

#include <map>
#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class MaterialPropertys;

class MaterialProperty {
public:
  MaterialProperty(MaterialPropertys& matPropertys);
  
  ~MaterialProperty();
  
  void load(const YAML::Node & node);
  
  virtual void breadboard(){}

  Simulation *root();
  EquationSystems *parent();

  MaterialPropertys &matPropertys_;
};


} // namespace nalu
} // namespace Sierra

#endif
