/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TemperaturePropAlgorithm_h
#define TemperaturePropAlgorithm_h

#include <Algorithm.h>

// standard c++
#include <string>

namespace stk {
namespace mesh {
class FieldBase;
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;
class PropertyEvaluator;

class TemperaturePropAlgorithm : public Algorithm
{
public:

  TemperaturePropAlgorithm(
    Realm & realm,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * prop,
    PropertyEvaluator *propEvaluator,
    std::string tempName = "temperature");

  virtual ~TemperaturePropAlgorithm() {}

  virtual void execute();

  stk::mesh::FieldBase *prop_;
  PropertyEvaluator *propEvaluator_;
  stk::mesh::FieldBase *temperature_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
